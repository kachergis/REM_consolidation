#!/usr/bin/env bash
# Publish a code-only copy of this repo's history to the public REM_consolidation repo.
#
# This repo (and its private origin) also holds reviews, slides and paper PDFs. This script
# clones the current branch to a temp dir, keeps only the ALLOW paths below in every commit,
# verifies that nothing else remains, and shows how the result relates to the public branch.
# With --push it force-updates that branch, which is a generated mirror: never commit to it
# directly. Filtering preserves dates, so unchanged commits keep the same hashes across runs.
#
# Usage: scripts/publish_public.sh          # dry run: build and verify only
#        scripts/publish_public.sh --push   # build, verify, push
set -euo pipefail

PUBLIC_URL=https://github.com/kachergis/REM_consolidation.git
PUBLIC_BRANCH=reimplementation
# Extended regex over repo-relative paths; anything not matching is removed from history.
ALLOW='^(rem_sleep/|modeling/|scripts/|\.gitignore$|cogsci\.sty$|holo_sleep_refs(_abbreviated)?\.bib$|androids_dream_[^/]*\.tex$|final ICCM/[^/]*\.tex$|final_graph_dist1_50reps\.pdf$)'

repo=$(git rev-parse --show-toplevel)
branch=$(git -C "$repo" branch --show-current)
if [ -n "$(git -C "$repo" status --porcelain)" ]; then
  echo "Working tree is not clean; commit or stash first." >&2
  exit 1
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
git clone -q --single-branch --branch "$branch" "$repo" "$tmp/mirror"
cd "$tmp/mirror"
git config core.quotepath off

export ALLOW
FILTER_BRANCH_SQUELCH_WARNING=1 git filter-branch --index-filter '
  git -c core.quotepath=off ls-files | grep -Ev "$ALLOW" |
    while IFS= read -r f; do git rm -q --cached -- "$f"; done
' HEAD >/dev/null 2>&1

leaked=$(git log HEAD --name-only --format= | sort -u | grep -v '^$' | grep -Ev "$ALLOW" || true)
if [ -n "$leaked" ]; then
  printf 'Refusing: paths outside the allowlist remain in history:\n%s\n' "$leaked" >&2
  exit 1
fi

echo "Filtered history of '$branch' ($(git rev-list --count HEAD) commits; all paths allowlisted):"
git log --oneline HEAD
echo
echo "Non-code paths that will be public:"
git log HEAD --name-only --format= | sort -u | grep -Ev '^(rem_sleep|modeling|scripts)/' | sed 's/^/  /'

remote_sha=$(git ls-remote "$PUBLIC_URL" "refs/heads/$PUBLIC_BRANCH" | cut -f1)
echo
if [ -z "$remote_sha" ]; then
  echo "Public branch '$PUBLIC_BRANCH' does not exist yet; it would be created."
elif [ "$remote_sha" = "$(git rev-parse HEAD)" ]; then
  echo "Public branch '$PUBLIC_BRANCH' is already up to date."
  exit 0
else
  git fetch -q "$PUBLIC_URL" "$PUBLIC_BRANCH"
  if git merge-base --is-ancestor FETCH_HEAD HEAD; then
    echo "Public branch is behind; this is a fast-forward:"
  else
    echo "Public branch has diverged; it will be force-updated. Commits it would lose:"
    git log --oneline HEAD..FETCH_HEAD | sed 's/^/  /'
    echo "New commits:"
  fi
  git log --oneline FETCH_HEAD..HEAD | sed 's/^/  /'
fi

if [ "${1:-}" != "--push" ]; then
  echo
  echo "Dry run; rerun with --push to publish."
  exit 0
fi
git push --force-with-lease="refs/heads/$PUBLIC_BRANCH:$remote_sha" "$PUBLIC_URL" \
  "HEAD:refs/heads/$PUBLIC_BRANCH"
