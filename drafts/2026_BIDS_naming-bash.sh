#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# bids_utils.bash  —  BIDS filename construction and entity injection
#
# Source: BIDS specification v1.11.1
#   Raw entities:         https://bids-specification.readthedocs.io/en/stable/appendices/entity-table.html
#   Derivative entities:  https://bids-specification.readthedocs.io/en/stable/derivatives/imaging.html
#
# Requires: bash 4.3+  (local -n nameref, [[ -v ]] key-existence test)
# ──────────────────────────────────────────────────────────────────────────────

# ── Canonical entity order ────────────────────────────────────────────────────
# Raw data entities (BIDS v1.11.1 entity table), followed by derivative-only
# entities.  desc is always last before the suffix, per spec.
#
# Entity      Full name                        Notes
# ----------  -------------------------------  --------------------------------
# sub         Subject                          all modalities
# ses         Session                          all modalities
# sample      Sample                           micr
# hemi        Hemisphere (L|R)                 surface files
# task        Task                             func, eeg, ieeg, meg, beh…
# tracksys    Tracking system                  motion
# acq         Acquisition                      all modalities
# ce          Contrast enhancing agent         anat, fmap
# trc         Tracer                           pet
# stain       Stain                            micr
# nuc         Nucleus                          mrs
# voi         Volume of interest               mrs
# rec         Reconstruction                   anat, dwi, fmap, perf, pet, mrs
# dir         Phase-encoding direction         dwi, func, fmap, perf
# run         Run index                        all modalities
# mod         Corresponding modality           anat (defacemask), noRF
# echo        Echo index                       anat, func, fmap, mrs
# flip        Flip angle index                 anat, fmap
# inv         Inversion time index             anat, fmap, mrs
# mt          Magnetization transfer           anat, fmap
# part        Part (mag/phase/real/imag)       anat, dwi, func, fmap
# proc        Processed on device              meg
# space       Coordinate space                 eeg, ieeg, meg + all derivs
# split       Split index                      meg
# recording   Recording label                  func, perf, dwi, anat…
# chunk       Chunk index                      anat, dwi, micr
# res         Resolution                       derivatives
# den         Surface density                  derivatives
# from        Source space (transforms)        derivatives
# to          Target space (transforms)        derivatives
# label       Anatomical structure label       derivatives (masks)
# seg         Segmentation label               derivatives
# atlas       Atlas label                      derivatives
# desc        Description (free-form)          derivatives — always last
# ──────────────────────────────────────────────────────────────────────────────
_BIDS_ORDER=(
    sub ses sample hemi task tracksys
    acq ce trc stain nuc voi
    rec dir run mod
    echo flip inv mt part
    proc space split recording chunk
    res den from to
    label seg atlas desc
)

# ──────────────────────────────────────────────────────────────────────────────
# _bids_rebuild <kv_var> <suffix> <ext> <out_var>
#
# Internal helper — pure bash, no I/O, no subshell, no fork.
# Iterates _BIDS_ORDER, pulls matching entries from the associative array
# <kv_var>, appends <suffix>, joins with underscores, appends <ext>, and
# stores the result in <out_var>.
#
# Double-underscore variable names prevent nameref circular-reference errors
# when callers happen to use the same local names (kv, suffix, ext, result…).
# ──────────────────────────────────────────────────────────────────────────────
_bids_rebuild() {
    local -n __kv="$1"
    local     __suffix="$2"
    local     __ext="$3"
    local -n  __result="$4"
    local -a  __parts
    local     __key
    for __key in "${_BIDS_ORDER[@]}"; do
        [[ -v __kv[$__key] ]] && __parts+=("${__kv[$__key]}")
    done
    [[ -n "$__suffix" ]] && __parts+=("$__suffix")
    local IFS=_
    __result="${__parts[*]}${__ext}"
}

# ──────────────────────────────────────────────────────────────────────────────
# bids_name  --<entity> <value> ...  --suffix <suffix>  [--ext <ext>]
#
# Build a BIDS filename from scratch using entity flags.
# Any key from _BIDS_ORDER is accepted as --<key>.
# --suffix is required; --ext is optional.
#
# Examples:
#   bids_name --sub 01 --ses 01 --task rest --run 1 \
#             --suffix bold --ext .nii.gz
#   → sub-01_ses-01_task-rest_run-1_bold.nii.gz
#
#   bids_name --sub 01 --hemi L --space fsLR --den 32k \
#             --suffix pial --ext .surf.gii
#   → sub-01_hemi-L_space-fsLR_den-32k_pial.surf.gii
#
#   bids_name --sub 01 --from T1w --to MNI152NLin6Asym \
#             --suffix xfm --ext .mat
#   → sub-01_from-T1w_to-MNI152NLin6Asym_xfm.mat
# ──────────────────────────────────────────────────────────────────────────────
bids_name() {
    local -A kv
    local suffix="" ext="" key result

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --suffix) suffix="$2";                           shift 2 ;;
            --ext)    ext="$2";                              shift 2 ;;
            --?*)     key="${1#--}"; kv[$key]="${key}-${2}"; shift 2 ;;
            *) printf 'bids_name: unknown argument: %s\n' "$1" >&2; return 1 ;;
        esac
    done

    [[ -z "$suffix" ]] && {
        printf 'bids_name: --suffix is required\n' >&2; return 1
    }

    _bids_rebuild kv "$suffix" "$ext" result
    printf '%s\n' "$result"
}

# ──────────────────────────────────────────────────────────────────────────────
# bids_add_entity  <file>  [--<entity> <value>]...
#
# Add or replace entities in an existing BIDS filename, reordering all
# entities into canonical BIDS order.
#
# <file> may be a bare filename or a full path; the directory prefix is
# preserved unchanged in the output.
#
# Extension is stripped at the first dot of the filename component, so
# multi-part extensions are all handled correctly:
#   .nii.gz   .surf.gii   .label.gii   .func.gii   .omat   .mat   …
#
# Optional overrides:
#   --suffix <label>   replace the existing suffix
#   --ext    <ext>     replace the extension (include the leading dot)
#
# Examples:
#   bids_add_entity sub-01_task-rest_bold.nii.gz --rec norm
#   → sub-01_task-rest_rec-norm_bold.nii.gz
#
#   bids_add_entity sub-01_echo-1_ses-01_bold.nii.gz --desc preproc
#   → sub-01_ses-01_echo-1_desc-preproc_bold.nii.gz   (reordered)
#
#   bids_add_entity sub-01_T1w.surf.gii --space MNI152 --desc preproc
#   → sub-01_space-MNI152_desc-preproc_T1w.surf.gii
#
#   bids_add_entity sub-01_from-T1w_to-MNI_xfm.omat --desc affine
#   → sub-01_from-T1w_to-MNI_desc-affine_xfm.omat
#
#   bids_add_entity /data/sub-01/func/sub-01_task-rest_bold.nii.gz \
#                   --rec norm --run 2
#   → /data/sub-01/func/sub-01_task-rest_rec-norm_run-2_bold.nii.gz
# ──────────────────────────────────────────────────────────────────────────────
bids_add_entity() {
    local f="$1"; shift

    # ── Separate directory prefix (pure string ops, no subshell) ──────────────
    local dir="" fname="$f"
    [[ "$f" == */* ]] && { dir="${f%/*}/"; fname="${f##*/}"; }

    # ── Strip extension at the first dot of the filename component ────────────
    local base ext
    if [[ "$fname" == *.* ]]; then
        base="${fname%%.*}"
        ext=".${fname#*.}"
    else
        base="$fname"
        ext=""
    fi

    # ── Parse existing entities and suffix ────────────────────────────────────
    local -A kv
    local -a parts
    local suffix="" p key
    IFS='_' read -ra parts <<< "$base"
    for p in "${parts[@]}"; do
        if [[ "$p" == *-* ]]; then
            key="${p%%-*}"; kv[$key]="$p"
        else
            suffix="$p"     # bare token (no hyphen) = BIDS suffix
        fi
    done

    # ── Apply new / overriding entities ───────────────────────────────────────
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --suffix) suffix="$2";                           shift 2 ;;
            --ext)    ext="$2";                              shift 2 ;;
            --?*)     key="${1#--}"; kv[$key]="${key}-${2}"; shift 2 ;;
            *) printf 'bids_add_entity: unknown argument: %s\n' "$1" >&2; return 1 ;;
        esac
    done

    local result
    _bids_rebuild kv "$suffix" "$ext" result
    printf '%s\n' "${dir}${result}"
}