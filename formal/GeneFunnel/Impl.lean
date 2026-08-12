import GeneFunnel.Spec

namespace GeneFunnel.Impl

def allNonnegative : List Rat → Bool
  | [] => true
  | value :: values => decide (0 ≤ value) && allNonnegative values

def below (mean value : Rat) : Bool := value < mean

def optimizedScore (values : List Rat) : Option Rat :=
  if values.length < 2 then none else
    let size : Rat := values.length
    let total := values.sum
    let mean := total / size
    let lows := values.filter (below mean)
    some (((size - 1 - lows.length) * total + size * lows.sum) /
      (size - 1))

def run : GeneFunnel.Spec.Program := fun cell =>
  let values := GeneFunnel.Spec.observed cell
  if allNonnegative values then
    .ok (optimizedScore values)
  else
    .error .negativeValue

end GeneFunnel.Impl
