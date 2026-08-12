import Std

namespace GeneFunnel.Spec

inductive Error where
  | negativeValue
  deriving DecidableEq, Repr

abbrev Cell := List (Option Rat)
abbrev Program := Cell → Except Error (Option Rat)

/-- Missing entries are absent; observed zeroes remain in the list. -/
def observed (cell : Cell) : List Rat := cell.filterMap id

def Valid (cell : Cell) : Prop :=
  ∀ value ∈ observed cell, 0 ≤ value

instance validDecidable (cell : Cell) : Decidable (Valid cell) := by
  unfold Valid
  exact List.decidableBAll (fun value : Rat => 0 ≤ value) (observed cell)

/-- Exact GeneFunnel equation; fewer than two observations have no score. -/
def directScore (values : List Rat) : Option Rat :=
  if values.length < 2 then none else
    let size : Rat := values.length
    let total := values.sum
    let mean := total / size
    let deviation := (values.map fun value => (value - mean).abs).sum
    some (total - size / (2 * (size - 1)) * deviation)

/-- Success and rejection are both sound and complete. -/
def Contract (program : Program) : Prop :=
  (∀ cell, Valid cell →
    program cell = .ok (directScore (observed cell))) ∧
  (∀ cell, ¬Valid cell → program cell = .error .negativeValue)

def BoundsContract (program : Program) : Prop :=
  ∀ cell output, program cell = .ok (some output) →
    0 ≤ output ∧ output ≤ (observed cell).sum

def reference : Program := fun cell =>
  if Valid cell then
    .ok (directScore (observed cell))
  else
    .error .negativeValue

theorem contract_inhabited : ∃ program, Contract program := by
  refine ⟨reference, ?_⟩
  constructor
  · intro cell h
    simp [reference, h]
  · intro cell h
    simp [reference, h]

end GeneFunnel.Spec
