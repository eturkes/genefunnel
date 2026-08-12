#!/bin/bash -p

if [[ $- != *p* ]] ||
    /usr/bin/env | /usr/bin/grep -q '^BASH_FUNC_'; then
  clean_path=${PATH:-/usr/local/bin:/usr/bin:/bin}
  /usr/bin/env -i HOME="${HOME:?HOME is required}" PATH="$clean_path" \
    /bin/bash --noprofile --norc -p "$0" "$@"
else

set -Eeuo pipefail
IFS=$' \t\n'
umask 022
unset -f -- "$(compgen -A function)" 2>/dev/null || true
for variable_name in BASH_ENV ENV ELAN_TOOLCHAIN LD_AUDIT LD_LIBRARY_PATH \
    LD_PRELOAD; do
  unset "$variable_name"
done
for variable_name in ${!DYLD_@} ${!LAKE_@} ${!LEAN_@}; do
  unset "$variable_name"
done

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly script_dir
project_root="$(cd -- "$script_dir/.." && pwd -P)"
readonly project_root
readonly input_manifest="formal/verification-inputs.txt"
readonly digest_manifest="formal/SHA256SUMS"
readonly expected_toolchain_hash="2bdc48adfa58d0017e538a0ad117c5d73d35deec879978f909406a80c8037273"
readonly expected_lakefile_hash="16e11432532227761f4a2ae3c09d9cdbdade1e13743ed7b0cbf3ee93d91f3eec"
readonly expected_lake_manifest_hash="fad6644ae97d354015b2d1ebe6ce24bc2d979ad2769a8067ffca46a3376e82df"
readonly expected_lean_commit="f3b06c705e6c85f5314019d5d3baab0fec5b580c"
readonly command_timeout=600
readonly output_limit=8388608
readonly output_blocks=$((output_limit / 1024))
readonly expected_inputs=(
  formal/GeneFunnel.lean
  formal/GeneFunnel/Audit.lean
  formal/GeneFunnel/Bounds.lean
  formal/GeneFunnel/Impl.lean
  formal/GeneFunnel/Proof.lean
  formal/GeneFunnel/Spec.lean
  formal/verification-inputs.txt
  formal/verify.sh
  inst/SCIENTIFIC_SPEC.md
  lake-manifest.json
  lakefile.toml
  lean-toolchain
  src/calculateScores.cpp
)
readonly expected_sources=(
  formal/GeneFunnel.lean
  formal/GeneFunnel/Audit.lean
  formal/GeneFunnel/Bounds.lean
  formal/GeneFunnel/Impl.lean
  formal/GeneFunnel/Proof.lean
  formal/GeneFunnel/Spec.lean
)
readonly model_sources=(
  formal/GeneFunnel.lean
  formal/GeneFunnel/Bounds.lean
  formal/GeneFunnel/Impl.lean
  formal/GeneFunnel/Proof.lean
  formal/GeneFunnel/Spec.lean
)
readonly audit_forbidden="(?<![A-Za-z0-9_'])(admit|axiom|constant|csimp|extern|implemented_by|opaque|partial_fixpoint|skipKernelTC|sorry|unsafe)(?![A-Za-z0-9_'])"
readonly model_forbidden="(?<![A-Za-z0-9_'])(admit|axiom|constant|csimp|extern|implemented_by|opaque|partial|partial_fixpoint|skipKernelTC|sorry|unsafe)(?![A-Za-z0-9_'])"
cd -- "$project_root"

fail() {
  printf 'formal verification: %s\n' "$*" >&2
  exit 1
}

for command_name in awk chmod cmp cp dirname head lake mkdir mktemp mv \
    realpath rg rm sha256sum sort stat timeout wc; do
  command_path="$(type -P -- "$command_name" || true)"
  [[ -n $command_path && -x $command_path ]] ||
    fail "required executable unavailable: $command_name"
done
cmd_awk="$(type -P awk)"
cmd_chmod="$(type -P chmod)"
cmd_cmp="$(type -P cmp)"
cmd_cp="$(type -P cp)"
cmd_dirname="$(type -P dirname)"
cmd_head="$(type -P head)"
cmd_lake="$(type -P lake)"
cmd_mkdir="$(type -P mkdir)"
cmd_mktemp="$(type -P mktemp)"
cmd_mv="$(type -P mv)"
cmd_realpath="$(type -P realpath)"
cmd_rg="$(type -P rg)"
cmd_rm="$(type -P rm)"
cmd_sha256sum="$(type -P sha256sum)"
cmd_sort="$(type -P sort)"
cmd_stat="$(type -P stat)"
cmd_timeout="$(type -P timeout)"
cmd_wc="$(type -P wc)"
readonly cmd_awk cmd_chmod cmd_cmp cmd_cp cmd_dirname cmd_head cmd_lake
readonly cmd_mkdir cmd_mktemp cmd_mv cmd_realpath cmd_rg cmd_rm cmd_sha256sum
readonly cmd_sort cmd_stat cmd_timeout cmd_wc

validate_inputs() {
  local base=$1 path resolved
  local -A seen=()
  cd -- "$base"
  [[ -f $input_manifest && ! -L $input_manifest ]] ||
    fail "missing regular input manifest"
  mapfile -t inputs < "$input_manifest"
  [[ $(printf '%s\n' "${inputs[@]}") == \
    "$(printf '%s\n' "${expected_inputs[@]}")" ]] ||
    fail "verification input inventory changed; update the reviewed boundary"
  for path in "${inputs[@]}"; do
    [[ $path =~ ^[A-Za-z0-9_][A-Za-z0-9._/-]*$ && \
      $path != /* && $path != *//* && $path != */./* && \
      $path != */. && $path != ../* && $path != */../* && $path != */.. ]] ||
      fail "unsafe input path: $path"
    [[ -f $path && ! -L $path ]] || fail "missing regular input: $path"
    resolved="$($cmd_realpath -e -- "$path")"
    [[ $resolved == "$base/$path" ]] || fail "input path is not canonical: $path"
    [[ -z ${seen[$path]+present} ]] || fail "duplicate input: $path"
    seen[$path]=1
  done
}

source_inventory() {
  local base=$1 output sorted
  cd -- "$base"
  output="$(LC_ALL=C "$cmd_rg" --no-config --files -uu formal -g '*.lean')" ||
    fail "cannot enumerate Lean sources"
  sorted="$(LC_ALL=C "$cmd_sort" <<<"$output")" ||
    fail "cannot sort Lean source inventory"
  [[ $sorted == "$(printf '%s\n' "${expected_sources[@]}")" ]] ||
    fail "Lean source inventory changed; update the gate's module boundary"
  [[ ! -e lakefile.lean && ! -L lakefile.lean ]] ||
    fail "lakefile.lean can bypass the reviewed declarative configuration"
  [[ ! -e .lake/package-overrides.json && \
    ! -L .lake/package-overrides.json ]] || fail "Lake overrides are not accepted"
}

hash_inputs() {
  local base=$1 output=$2
  cd -- "$base"
  "$cmd_sha256sum" -- "${expected_inputs[@]}" > "$output" ||
    fail "cannot hash verification inputs"
}

check_fixed_hash() {
  local path=$1 expected=$2 observed
  observed="$($cmd_sha256sum -- "$path")" || fail "cannot hash $path"
  [[ ${observed%% *} == "$expected" ]] ||
    fail "$path differs from the reviewed build closure"
}

validate_snapshot() {
  local base=$1 expected_manifest pattern scan_rc reviewed_lines
  validate_inputs "$base"
  source_inventory "$base"
  cd -- "$base"
  expected_manifest="$base/expected-inputs.txt"
  printf '%s\n' "${expected_inputs[@]}" > "$expected_manifest"
  "$cmd_cmp" --silent -- "$input_manifest" "$expected_manifest" ||
    fail "input manifest must exactly encode the reviewed inventory"
  check_fixed_hash lean-toolchain "$expected_toolchain_hash"
  check_fixed_hash lakefile.toml "$expected_lakefile_hash"
  check_fixed_hash lake-manifest.json "$expected_lake_manifest_hash"
  for pattern in model audit; do
    set +e
    if [[ $pattern == model ]]; then
      LC_ALL=C "$cmd_rg" --no-config --pcre2 -n "$model_forbidden" \
        "${model_sources[@]}"
    else
      LC_ALL=C "$cmd_rg" --no-config --pcre2 -n "$audit_forbidden" \
        formal/GeneFunnel/Audit.lean
    fi
    scan_rc=$?
    set -e
    case $scan_rc in
      0) fail "forbidden Lean declaration spelling found in $pattern sources" ;;
      1) ;;
      *) fail "$pattern source-policy scan failed with exit $scan_rc" ;;
    esac
  done
  reviewed_lines="$($cmd_awk 'NF { count++ } END { print count + 0 }' \
    formal/GeneFunnel/Spec.lean)" || fail "cannot measure the trusted Lean surface"
  ((reviewed_lines <= 80)) || fail "trusted Lean surface exceeds 80 nonblank lines"
  printf '%s\n' "$reviewed_lines"
}

run_bounded() {
  local label=$1 log rc bytes
  shift
  log="$snapshot/${label//[^A-Za-z0-9_.-]/_}.log"
  set +e
  (
    ulimit -f "$output_blocks"
    "$cmd_timeout" --kill-after=5s "${command_timeout}s" "$@"
  ) > "$log" 2>&1
  rc=$?
  set -e
  bytes="$($cmd_wc -c < "$log")"
  if ((bytes >= output_limit && rc != 0)); then
    "$cmd_head" -c "$output_limit" -- "$log" >&2
    fail "$label exceeded the ${output_limit}-byte output limit"
  fi
  if ((rc != 0)); then
    "$cmd_head" -c "$output_limit" -- "$log" >&2
    fail "$label failed with exit $rc"
  fi
  "$cmd_head" -c "$output_limit" -- "$log"
}

validate_inputs "$project_root"
source_inventory "$project_root"

snapshot="$($cmd_mktemp -d -t genefunnel-formal.XXXXXXXX)"
readonly snapshot
published_digest=""
# shellcheck disable=SC2329
cleanup() {
  local rc=$?
  trap - EXIT
  if [[ -n ${published_digest:-} && -f $published_digest ]]; then
    "$cmd_rm" -f -- "$published_digest"
  fi
  if [[ -n ${snapshot:-} && -d $snapshot ]]; then
    "$cmd_rm" -rf -- "$snapshot"
  fi
  exit "$rc"
}
trap cleanup EXIT

if [[ ${1:-} != --update-hashes ]]; then
  (($# == 0)) || fail "usage: formal/verify.sh [--update-hashes]"
  [[ -f $digest_manifest && ! -L $digest_manifest ]] ||
    fail "missing regular digest record: $digest_manifest"
  [[ $($cmd_realpath -e -- "$digest_manifest") == \
    "$project_root/$digest_manifest" ]] || fail "digest record path is not canonical"
  [[ $($cmd_stat -c '%a' -- "$digest_manifest") == 644 ]] ||
    fail "digest record must have mode 0644"
  "$cmd_cp" --preserve=mode,timestamps -- "$digest_manifest" \
    "$snapshot/initial-SHA256SUMS"
  hash_inputs "$project_root" "$snapshot/root-before-SHA256SUMS"
  "$cmd_cmp" --silent -- "$snapshot/initial-SHA256SUMS" \
    "$snapshot/root-before-SHA256SUMS" ||
    fail "inputs differ from the complete digest record"
else
  (($# == 1)) || fail "usage: formal/verify.sh [--update-hashes]"
fi

for path in "${expected_inputs[@]}"; do
  "$cmd_mkdir" -p -- "$snapshot/$($cmd_dirname -- "$path")"
  "$cmd_cp" --preserve=mode,timestamps -- "$project_root/$path" "$snapshot/$path"
done
hash_inputs "$snapshot" "$snapshot/snapshot-SHA256SUMS"
if [[ ${1:-} != --update-hashes ]]; then
  "$cmd_cmp" --silent -- "$snapshot/initial-SHA256SUMS" \
    "$snapshot/snapshot-SHA256SUMS" || fail "snapshot differs from the digest record"
fi
reviewed_lines="$(validate_snapshot "$snapshot")"
readonly reviewed_lines

if [[ ${1:-} == --update-hashes ]]; then
  cd -- "$project_root"
  validate_inputs "$project_root"
  source_inventory "$project_root"
  hash_inputs "$project_root" "$snapshot/root-after-SHA256SUMS"
  "$cmd_cmp" --silent -- "$snapshot/snapshot-SHA256SUMS" \
    "$snapshot/root-after-SHA256SUMS" || fail "inputs changed while updating hashes"
  published_digest="$($cmd_mktemp -p formal .SHA256SUMS.XXXXXXXX)"
  "$cmd_cp" -- "$snapshot/snapshot-SHA256SUMS" "$published_digest"
  "$cmd_chmod" 0644 "$published_digest"
  hash_inputs "$project_root" "$snapshot/root-final-SHA256SUMS"
  "$cmd_cmp" --silent -- "$snapshot/snapshot-SHA256SUMS" \
    "$snapshot/root-final-SHA256SUMS" || fail "inputs changed before publishing hashes"
  "$cmd_mv" -f -- "$published_digest" "$digest_manifest"
  published_digest=""
  printf 'formal verification: updated %s\n' "$digest_manifest"
  exit 0
fi

cd -- "$snapshot"
[[ ! -e .lake && ! -L .lake ]] || fail "snapshot unexpectedly contains Lake state"
lean_version="$(run_bounded toolchain "$cmd_lake" env lean --version)"
readonly lean_version
readonly lean_version_pattern="^Lean \\(version 4\\.32\\.2, [^,]+, commit ${expected_lean_commit}, Release\\)$"
[[ $lean_version =~ $lean_version_pattern ]] ||
  fail "unexpected Lean build: $lean_version"
run_bounded build "$cmd_lake" --rehash --reconfigure --no-cache build --wfail \
  +GeneFunnel +GeneFunnel.Proof +GeneFunnel.Bounds
run_bounded audit "$cmd_lake" env lean formal/GeneFunnel/Audit.lean
run_bounded proof-replay "$cmd_lake" env leanchecker --fresh GeneFunnel.Proof
run_bounded bounds-replay "$cmd_lake" env leanchecker --fresh GeneFunnel.Bounds
hash_inputs "$snapshot" "$snapshot/snapshot-final-SHA256SUMS"
"$cmd_cmp" --silent -- "$snapshot/initial-SHA256SUMS" \
  "$snapshot/snapshot-final-SHA256SUMS" ||
  fail "verification inputs changed inside the snapshot"

cd -- "$project_root"
validate_inputs "$project_root"
source_inventory "$project_root"
"$cmd_cmp" --silent -- "$digest_manifest" "$snapshot/initial-SHA256SUMS" ||
  fail "digest record changed during verification"
hash_inputs "$project_root" "$snapshot/root-final-SHA256SUMS"
"$cmd_cmp" --silent -- "$snapshot/initial-SHA256SUMS" \
  "$snapshot/root-final-SHA256SUMS" || fail "inputs changed during verification"
printf 'formal verification: passed (%s trusted lines; two obligations)\n' \
  "$reviewed_lines"
fi
