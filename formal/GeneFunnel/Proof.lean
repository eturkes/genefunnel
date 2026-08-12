import GeneFunnel.Impl
import Init.Grind
import Init.Omega

namespace GeneFunnel.Proof

open GeneFunnel

def low (mean : Rat) (values : List Rat) :=
  values.filter (Impl.below mean)

def high (mean : Rat) (values : List Rat) :=
  values.filter fun value => !Impl.below mean value

theorem partition_sum (mean : Rat) (values : List Rat) :
    (low mean values).sum + (high mean values).sum = values.sum := by
  induction values with
  | nil =>
      simp only [low, high, List.filter_nil, List.sum_nil, Rat.zero_add]
  | cons value values ih =>
      simp only [low, high, List.filter_cons, Impl.below] at ih ⊢
      split <;> grind

theorem partition_length (mean : Rat) (values : List Rat) :
    (low mean values).length + (high mean values).length = values.length := by
  induction values with
  | nil => rfl
  | cons value values ih =>
      simp only [low, high, List.filter_cons, Impl.below] at ih ⊢
      by_cases h : value < mean <;> simp [h] at * <;> omega

theorem sum_abs (mean : Rat) (values : List Rat) :
    (values.map fun value => (value - mean).abs).sum =
      ((low mean values).length : Rat) * mean - (low mean values).sum +
      ((high mean values).sum - (high mean values).length * mean) := by
  induction values with
  | nil =>
      simp only [low, high, List.filter_nil, List.map_nil, List.sum_nil,
        List.length_nil]
      grind
  | cons value values ih =>
      simp only [List.map_cons, List.sum_cons, low, high, List.filter_cons,
        Impl.below] at ih ⊢
      by_cases h : value < mean
      · have habs : (value - mean).abs = mean - value := by
          rw [Rat.abs_of_nonpos (by grind)]
          grind
        simp only [h, decide_true, ↓reduceIte, Bool.not_true,
          Bool.false_eq_true, List.sum_cons, List.length_cons]
        rw [habs, ih]
        grind
      · have habs : (value - mean).abs = value - mean := by
          rw [Rat.abs_of_nonneg (by grind)]
        simp only [h, decide_false, ↓reduceIte, Bool.not_false,
          List.sum_cons, List.length_cons]
        rw [habs, ih]
        grind

theorem score_identity (values : List Rat) :
    Impl.optimizedScore values = Spec.directScore values := by
  unfold Impl.optimizedScore Spec.directScore
  split
  · rfl
  · rename_i hsize
    simp only
    let size : Rat := values.length
    let total := values.sum
    let mean := total / size
    let lows := low mean values
    let highs := high mean values
    have hsizeNat : 2 ≤ values.length := by omega
    have hsizePos : 0 < size := by
      simp [size, Rat.natCast_pos]
      omega
    have hdenomPos : 0 < size - 1 := by
      have hcast : (1 : Rat) < size := by
        rw [show (1 : Rat) = (1 : Nat) by rfl]
        exact Rat.natCast_lt_natCast.mpr (by omega)
      grind
    have hmean : mean * size = total :=
      Rat.div_mul_cancel (Rat.ne_of_gt hsizePos)
    have hparts : lows.sum + highs.sum = total :=
      partition_sum mean values
    have hlensNat : lows.length + highs.length = values.length :=
      partition_length mean values
    have hlens : (lows.length : Rat) + highs.length = size := by
      simpa [size] using congrArg (fun count : Nat => (count : Rat)) hlensNat
    have habs := sum_abs mean values
    apply congrArg some
    change (((size - 1 - lows.length) * total + size * lows.sum) /
      (size - 1)) = total - size / (2 * (size - 1)) *
      (values.map fun value => (value - mean).abs).sum
    rw [habs]
    change (((size - 1 - lows.length) * total + size * lows.sum) /
      (size - 1)) = total - size / (2 * (size - 1)) *
      (lows.length * mean - lows.sum +
        (highs.sum - highs.length * mean))
    grind

theorem all_nonnegative (values : List Rat) :
    Impl.allNonnegative values = true ↔ ∀ value ∈ values, 0 ≤ value := by
  induction values with
  | nil => simp [Impl.allNonnegative]
  | cons value values ih => simp [Impl.allNonnegative, ih]

theorem implementation_correct : Spec.Contract Impl.run := by
  constructor
  · intro cell hvalid
    simp [Impl.run, (all_nonnegative _).mpr hvalid, score_identity]
  · intro cell hinvalid
    have hcheck : Impl.allNonnegative (Spec.observed cell) = false := by
      cases h : Impl.allNonnegative (Spec.observed cell)
      · rfl
      · exact False.elim (hinvalid ((all_nonnegative _).mp h))
    simp [Impl.run, hcheck]

end GeneFunnel.Proof
