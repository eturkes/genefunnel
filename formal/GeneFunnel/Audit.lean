import Lean.Compiler.ExternAttr
import Lean.Compiler.ImplementedByAttr
import Lean.Compiler.NoncomputableAttr
import Lean.Elab.Command
import Lean.Parser.Module
import Lean.Util.CollectAxioms

open Lean Elab Command

private def originOf? (env : Environment) (name : Name) : Option Name := do
  let index ← env.getModuleIdxFor? name
  env.header.moduleNames[index]?

private def assertOrigin
    (env : Environment) (name expected : Name) : CommandElabM Unit := do
  let some actual := originOf? env name
    | throwError m!"missing module origin for {name}"
  unless actual == expected do
    throwError m!"{name} originates in {actual}, expected {expected}"

private partial def auditClosure
    (env : Environment) (localModules trusted : Array Name) (trustedOnly : Bool)
    (pending seen : List Name) : CommandElabM Unit := do
  let some name := pending.head? | return
  let rest := pending.tail!
  if seen.contains name then
    auditClosure env localModules trusted trustedOnly rest seen
    return
  let seen := name :: seen
  let some origin := originOf? env name
    | auditClosure env localModules trusted trustedOnly rest seen
      return
  if !localModules.contains origin then
    auditClosure env localModules trusted trustedOnly rest seen
    return
  if trustedOnly && !trusted.contains origin then
    throwError m!"trusted semantics depend on {name} from {origin}"
  let some info := env.find? name
    | throwError m!"missing local declaration {name}"
  if (Lean.Compiler.getImplementedBy? env name).isSome then
    throwError m!"local declaration {name} has a runtime replacement"
  if Lean.isExtern env name then
    throwError m!"local declaration {name} has an external runtime implementation"
  if !trustedOnly && Lean.isNoncomputable env name then
    throwError m!"local executable dependency {name} is noncomputable"
  match info with
  | .axiomInfo _ => throwError m!"unsupported local declaration {name}"
  | .opaqueInfo _ => throwError m!"unsupported local declaration {name}"
  | .defnInfo value =>
      unless value.safety == .safe do
        throwError m!"local definition {name} is not total and safe"
  | .inductInfo value =>
      if value.isUnsafe then
        throwError m!"local inductive {name} is not safe"
  | .ctorInfo value =>
      if value.isUnsafe then
        throwError m!"local constructor {name} is not safe"
  | .recInfo value =>
      if value.isUnsafe then
        throwError m!"local recursor {name} is not safe"
  | _ => pure ()
  let pending := info.type.foldConsts rest fun dependency queue =>
    dependency :: queue
  let pending := match info.value? with
    | some value => value.foldConsts pending fun dependency queue =>
        dependency :: queue
    | none => pending
  auditClosure env localModules trusted trustedOnly pending seen

private def auditHeader
    (path : String) (localModules trustedModules : Array Name) :
    CommandElabM Unit := do
  let input ← IO.FS.readFile path
  let context := Parser.mkInputContext input path
  let (header, _, messages) ← Parser.parseHeader context
  if messages.hasErrors then
    throwError m!"cannot parse trusted header {path}"
  if let `(Parser.Module.header| $[module%$moduleTk?]? $[prelude]?
      $importsStx*) := header then
    for importStx in importsStx do
      if let `(Parser.Module.import| $[public%$pubTk?]? $[meta%$metaTk?]?
          import $[all%$allTk?]? $moduleName) := importStx then
        let imported := moduleName.getId
        if localModules.contains imported && !trustedModules.contains imported then
          throwError m!"trusted file {path} imports local module {imported}"
  else
    throwError m!"cannot decode trusted header {path}"

private def auditObligation
    (env : Environment) (localModules trustedModules : Array Name)
    (declaration declarationModule contract contractModule implementation
      implementationModule : Name) : CommandElabM Unit := do
  assertOrigin env declaration declarationModule
  assertOrigin env contract contractModule
  assertOrigin env implementation implementationModule
  let some declarationInfo := env.find? declaration
    | throwError m!"missing obligation {declaration}"
  unless declarationInfo.isTheorem do
    throwError m!"obligation {declaration} is not a theorem"
  let some contractInfo := env.find? contract
    | throwError m!"missing contract {contract}"
  unless contractInfo.isDefinition do
    throwError m!"contract {contract} is not a definition"
  let some implementationInfo := env.find? implementation
    | throwError m!"missing implementation {implementation}"
  unless implementationInfo.isDefinition do
    throwError m!"implementation {implementation} is not a definition"
  match declarationInfo.type.consumeMData with
  | .app (.const actualContract _) (.const actualImplementation _) =>
      unless actualContract == contract && actualImplementation == implementation do
        throwError m!"wrong obligation roots: {repr declarationInfo.type}"
  | _ =>
      throwError m!"obligation is not an exact unary contract application: \
        {repr declarationInfo.type}"
  auditClosure env localModules trustedModules true [contract] []
  auditClosure env localModules trustedModules false [implementation] []
  let previous ← getEnv
  setEnv env
  let observed ← Lean.collectAxioms declaration
  setEnv previous
  let expected : Array Name :=
    #[`Classical.choice, `Quot.sound, `propext]
  for name in observed.toList do
    unless expected.contains name do
      throwError m!"unexpected logical dependency {name} in {declaration}"

run_cmd do
  let env ← Lean.importModules #[{ module := `GeneFunnel.Bounds }] {}
  let localModules : Array Name := #[
    `GeneFunnel,
    `GeneFunnel.Audit,
    `GeneFunnel.Bounds,
    `GeneFunnel.Impl,
    `GeneFunnel.Proof,
    `GeneFunnel.Spec,
  ]
  let trustedModules : Array Name := #[`GeneFunnel.Spec]
  auditHeader "formal/GeneFunnel/Spec.lean" localModules trustedModules
  auditObligation env localModules trustedModules
    `GeneFunnel.Proof.implementation_correct `GeneFunnel.Proof
    `GeneFunnel.Spec.Contract `GeneFunnel.Spec
    `GeneFunnel.Impl.run `GeneFunnel.Impl
  auditObligation env localModules trustedModules
    `GeneFunnel.Bounds.bounds_correct `GeneFunnel.Bounds
    `GeneFunnel.Spec.BoundsContract `GeneFunnel.Spec
    `GeneFunnel.Impl.run `GeneFunnel.Impl
