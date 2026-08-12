import GeneFunnel.Proof

namespace GeneFunnel.Bounds

open GeneFunnel

theorem sum_nonnegative {values : List Rat}
    (h : ∀ value ∈ values, 0 ≤ value) : 0 ≤ values.sum := by
  induction values with
  | nil => exact Rat.le_refl
  | cons value values ih =>
      exact Rat.add_nonneg (h value (by simp)) (ih (by
        intro other hother
        exact h other (by simp [hother])))

theorem low_sum_nonnegative {mean : Rat} {values : List Rat}
    (h : ∀ value ∈ values, 0 ≤ value) :
    0 ≤ (Proof.low mean values).sum := by
  apply sum_nonnegative
  intro value hvalue
  apply h value
  exact (List.mem_filter.mp hvalue).1

theorem sum_lt_mul_of_all_lt (mean : Rat) {values : List Rat}
    (hne : values ≠ []) (hall : ∀ value ∈ values, value < mean) :
    values.sum < (values.length : Rat) * mean := by
  induction values with
  | nil => contradiction
  | cons value values ih =>
      have hvalue : value < mean := hall value (by simp)
      by_cases htail : values = []
      · subst values
        simp only [List.sum_cons, List.sum_nil, List.length_cons,
          List.length_nil, Rat.natCast_add]
        grind
      · have hrest : values.sum < (values.length : Rat) * mean := by
          apply ih htail
          intro other hother
          exact hall other (by simp [hother])
        simp only [List.sum_cons, List.length_cons, Rat.natCast_add]
        grind

theorem exists_not_below_mean {values : List Rat} (hlen : 0 < values.length) :
    ∃ value ∈ values, ¬value < values.sum / (values.length : Rat) := by
  apply Decidable.byContradiction
  intro hnone
  have hall : ∀ value ∈ values,
      value < values.sum / (values.length : Rat) := by
    intro value hvalue
    by_cases hlt : value < values.sum / (values.length : Rat)
    · exact hlt
    · exact False.elim (hnone ⟨value, hvalue, hlt⟩)
  have hne : values ≠ [] := by
    intro heq
    subst values
    simp at hlen
  have hlt := sum_lt_mul_of_all_lt
    (values.sum / (values.length : Rat)) hne hall
  have hn : (values.length : Rat) ≠ 0 := by
    intro hz
    have : values.length = 0 := Rat.natCast_eq_zero_iff.mp hz
    omega
  have hmul : (values.length : Rat) *
      (values.sum / (values.length : Rat)) = values.sum := by
    rw [Rat.mul_comm]
    exact Rat.div_mul_cancel hn
  rw [hmul] at hlt
  exact Rat.lt_irrefl hlt

theorem low_length_lt {values : List Rat} (hlen : 0 < values.length) :
    (Proof.low (values.sum / (values.length : Rat)) values).length <
      values.length := by
  rw [Proof.low, List.length_filter_lt_length_iff_exists]
  obtain ⟨value, hvalue, hnot⟩ := exists_not_below_mean hlen
  exact ⟨value, hvalue, by simp [Impl.below, hnot]⟩

theorem low_count_le_previous {values : List Rat}
    (hlen : 2 ≤ values.length) :
    ((Proof.low (values.sum / (values.length : Rat)) values).length : Rat) ≤
      (values.length : Rat) - 1 := by
  have hlt := low_length_lt (values := values) (by omega)
  have hnat :
      (Proof.low (values.sum / (values.length : Rat)) values).length + 1 ≤
        values.length := by omega
  have hcast := Rat.natCast_le_natCast.mpr hnat
  simp only [Rat.natCast_add] at hcast
  grind

theorem optimized_output_nonnegative {values : List Rat} {output : Rat}
    (hvalues : ∀ value ∈ values, 0 ≤ value)
    (hscore : Impl.optimizedScore values = some output) : 0 ≤ output := by
  unfold Impl.optimizedScore at hscore
  split at hscore
  · contradiction
  · rename_i hsize
    simp only at hscore
    cases hscore
    let size : Rat := values.length
    let total := values.sum
    let mean := total / size
    let lows := Proof.low mean values
    have hlen : 2 ≤ values.length := by omega
    have hdenom : 0 < size - 1 := by
      have hcast : (1 : Rat) < size := by
        rw [show (1 : Rat) = (1 : Nat) by rfl]
        exact Rat.natCast_lt_natCast.mpr (by omega)
      grind
    have hcount : (lows.length : Rat) ≤ size - 1 :=
      low_count_le_previous hlen
    have hleft : 0 ≤ size - 1 - lows.length := by grind
    have htotal : 0 ≤ total := sum_nonnegative hvalues
    have hlow : 0 ≤ lows.sum := low_sum_nonnegative hvalues
    have hfirst := Rat.mul_nonneg hleft htotal
    have hsecond : 0 ≤ size * lows.sum :=
      Rat.mul_nonneg Rat.natCast_nonneg hlow
    have hnum := Rat.add_nonneg hfirst hsecond
    rw [Rat.div_def]
    exact Rat.mul_nonneg hnum (Rat.le_of_lt (Rat.inv_pos.mpr hdenom))

theorem deviation_nonnegative (mean : Rat) (values : List Rat) :
    0 ≤ (values.map fun value => (value - mean).abs).sum := by
  apply sum_nonnegative
  intro value hvalue
  obtain ⟨source, _, rfl⟩ := List.mem_map.mp hvalue
  exact Rat.abs_nonneg

theorem direct_output_le_sum {values : List Rat} {output : Rat}
    (hscore : Spec.directScore values = some output) : output ≤ values.sum := by
  unfold Spec.directScore at hscore
  split at hscore
  · contradiction
  · rename_i hsize
    simp only at hscore
    cases hscore
    let size : Rat := values.length
    let total := values.sum
    let mean := total / size
    have hlen : 2 ≤ values.length := by omega
    have hdenom : 0 < 2 * (size - 1) := by
      have hcast : (1 : Rat) < size := by
        rw [show (1 : Rat) = (1 : Nat) by rfl]
        exact Rat.natCast_lt_natCast.mpr (by omega)
      exact Rat.mul_pos (by decide) (by grind)
    have hcoef : 0 ≤ size / (2 * (size - 1)) := by
      rw [Rat.div_def]
      exact Rat.mul_nonneg Rat.natCast_nonneg
        (Rat.le_of_lt (Rat.inv_pos.mpr hdenom))
    have hdev := deviation_nonnegative mean values
    have hpenalty := Rat.mul_nonneg hcoef hdev
    grind

theorem bounds_correct : Spec.BoundsContract Impl.run := by
  intro cell output hrun
  unfold Impl.run at hrun
  dsimp only at hrun
  split at hrun
  · rename_i hcheck
    have hvalues := (Proof.all_nonnegative _).mp hcheck
    have hscore : Impl.optimizedScore (Spec.observed cell) = some output :=
      Except.ok.inj hrun
    constructor
    · exact optimized_output_nonnegative hvalues hscore
    · apply direct_output_le_sum
      rw [← Proof.score_identity]
      exact hscore
  · contradiction

end GeneFunnel.Bounds
