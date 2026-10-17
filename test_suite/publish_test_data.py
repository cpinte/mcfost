#!/usr/bin/env python3
"""Regenerate and publish the mcfost test-suite reference data.

Automates the workflow that was previously done by hand:

  1. Build mcfost single-threaded (mono=yes) for reproducible reference data.
  2. Run compute_test_suite.sh to (re)generate test_data/<model>/data_*.
  3. Package test_data/ into a dated, git-sha-tagged tarball, matching the
     existing naming convention (test_data_<sha>_<DD_MM_YYYY>.tgz).
  4. Upload the tarball to Google Drive via rclone, and update
     get_test_data.sh so the test suite fetches the new version.

Each stage can be skipped independently (e.g. to only re-package data you
already generated, or to test the publish step without rebuilding). By
default the script packages, uploads and rewrites get_test_data.sh; pass
--dry-run to print what the package/upload/rewrite steps would do instead,
since those touch the shared Google Drive artifact every CI run depends on.

Requirements for the upload step:
  - rclone (https://rclone.org/drive/), configured once with a remote that
    points at the shared Google Drive folder already hosting test_data_*.tgz:
        rclone config
    (name the remote anything; pass it via --remote NAME:path/to/folder,
    or set MCFOST_TEST_DATA_REMOTE in your environment.)
  - That folder must already be shared as "anyone with the link can view" so
    get_test_data.sh's anonymous download link keeps working for new files
    added to it (Drive normally applies a folder's sharing setting to new
    files uploaded into it, but double-check the first time).

Usage examples:
  # Full publish, dry-run first to see what would change:
  ./publish_test_data.py --remote gdrive:mcfost_test_data --dry-run

  # Actually do it, once you're happy with the dry-run output:
  ./publish_test_data.py --remote gdrive:mcfost_test_data

  # Package and upload the existing test_data/, without rebuilding or
  # re-running the calculations:
  ./publish_test_data.py --upload-only --remote ...

  # Upload a tarball that was already packaged (e.g. after a failed upload):
  ./publish_test_data.py --upload-only --tarball test_data_<sha>_<date>.tgz --remote ...

  # Only build + generate + tar, skip the upload entirely:
  ./publish_test_data.py --no-upload
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tarfile
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "src"
TEST_SUITE_DIR = REPO_ROOT / "test_suite"
TEST_DATA_DIR = TEST_SUITE_DIR / "test_data"
GET_TEST_DATA_SH = TEST_SUITE_DIR / "get_test_data.sh"


def run(cmd, cwd=None, env=None, dry_run=False):
    printable = " ".join(str(c) for c in cmd)
    print(f"$ {printable}" + (f"   (cwd={cwd})" if cwd else ""))
    if dry_run:
        print("  [dry-run: not executed]")
        return
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def git_sha():
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
    ).strip()


def build_mcfost(compiler, jobs, dry_run):
    print(f"\n== Building mcfost ({compiler}, single-threaded: mono=yes) ==")
    run(["make", f"{compiler}=yes", "mono=yes", f"-j{jobs}"], cwd=SRC_DIR, dry_run=dry_run)


def generate_reference_data(dry_run):
    print("\n== Regenerating reference data (compute_test_suite.sh) ==")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    run(["./compute_test_suite.sh"], cwd=TEST_SUITE_DIR, env=env, dry_run=dry_run)


def make_tarball(sha, tag_date, dry_run):
    name = f"test_data_{sha}_{tag_date}.tgz"
    path = TEST_SUITE_DIR / name
    print(f"\n== Packaging {name} ==")
    if dry_run:
        print(f"  [dry-run: would tar {TEST_DATA_DIR} -> {path}]")
        return path
    if not TEST_DATA_DIR.is_dir():
        sys.exit(f"error: {TEST_DATA_DIR} does not exist - run without --skip-generate first")
    with tarfile.open(path, "w:gz") as tar:
        tar.add(TEST_DATA_DIR, arcname="test_data")
    size_mb = path.stat().st_size / 1e6
    print(f"  wrote {path} ({size_mb:.1f} MB)")
    return path


def parse_drive_file_id(link):
    # rclone link returns either .../file/d/<id>/view or .../open?id=<id>
    match = re.search(r"/d/([\w-]+)", link) or re.search(r"[?&]id=([\w-]+)", link)
    return match.group(1) if match else None


def upload_rclone(tarball, remote, dry_run):
    print(f"\n== Uploading {tarball.name} to {remote} ==")
    if shutil.which("rclone") is None:
        print("  rclone not found - skipping upload.")
        print("  Install it and run `rclone config` once to set up a Google Drive remote:")
        print("  https://rclone.org/drive/")
        return None, None

    remote_path = f"{remote.rstrip('/')}/{tarball.name}"
    if dry_run:
        print(f"  [dry-run: would run `rclone copy {tarball} {remote}`]")
        print(f"  [dry-run: would run `rclone link {remote_path}` to get the shareable ID]")
        return None, None

    run(["rclone", "copy", str(tarball), remote], dry_run=dry_run)
    link = subprocess.check_output(["rclone", "link", remote_path], text=True).strip()
    file_id = parse_drive_file_id(link)
    if file_id is None:
        print(f"  Uploaded, but could not parse a file ID out of the link: {link}")
        print("  You'll need to grab the ID manually from Google Drive and update get_test_data.sh.")
    return file_id, link


def update_get_test_data_sh(new_name, new_id, dry_run):
    print(f"\n== Updating {GET_TEST_DATA_SH.name} ==")
    text = GET_TEST_DATA_SH.read_text()

    old_name_match = re.search(r'TEST_DATA_FILE="([^"]+)"', text)
    old_name = old_name_match.group(1) if old_name_match else None
    today = date.today().strftime("%d/%m/%Y")

    new_text = re.sub(r'TEST_DATA_FILE="[^"]+"', f'TEST_DATA_FILE="{new_name}"', text)
    new_text = re.sub(r'TEST_DATA_ID="[^"]+"', f'TEST_DATA_ID="{new_id}"', new_text)

    if old_name is not None:
        history_line = f"#{old_name} until {today}\n"
        marker = "#Previous versions:\n"
        if marker in new_text:
            new_text = new_text.replace(marker, marker + history_line, 1)
        else:
            new_text += f"\n\n{marker}{history_line}"

    if dry_run:
        print("  [dry-run: would rewrite get_test_data.sh as follows]")
        print("  ---")
        for line in new_text.splitlines():
            print("  " + line)
        print("  ---")
        return

    GET_TEST_DATA_SH.write_text(new_text)
    print(f"  updated TEST_DATA_FILE -> {new_name}")
    print(f"  updated TEST_DATA_ID   -> {new_id}")
    if old_name is not None:
        print(f"  recorded {old_name!r} in the Previous versions block")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--skip-build", action="store_true", help="Skip `make mono=yes`")
    ap.add_argument(
        "--skip-generate", action="store_true", help="Skip compute_test_suite.sh"
    )
    ap.add_argument("--no-upload", action="store_true", help="Only build+generate+tar")
    ap.add_argument(
        "--upload-only",
        action="store_true",
        help="Skip build and generate: package the existing test_data/ "
        "(or use --tarball) and upload it",
    )
    ap.add_argument(
        "--tarball",
        type=Path,
        help="Upload this existing tarball instead of packaging test_data/ "
        "(requires --upload-only)",
    )
    ap.add_argument(
        "--compiler", default="gfortran", choices=["gfortran", "ifort", "ifx"]
    )
    ap.add_argument("-j", "--jobs", default="1")
    ap.add_argument(
        "--remote",
        default=os.environ.get("MCFOST_TEST_DATA_REMOTE", ""),
        help="rclone remote:path for the shared Drive folder, "
        "e.g. gdrive:mcfost_test_data (or set MCFOST_TEST_DATA_REMOTE)",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what the package/upload/rewrite steps would do, without "
        "doing them (default: package, upload and rewrite get_test_data.sh)",
    )
    args = ap.parse_args()

    if args.tarball and not args.upload_only:
        ap.error("--tarball requires --upload-only")
    if args.upload_only and args.no_upload:
        ap.error("--upload-only and --no-upload are mutually exclusive")
    if args.upload_only:
        args.skip_build = args.skip_generate = True

    if args.dry_run:
        print("*** DRY RUN: no files will be uploaded or rewritten ***\n")

    if not args.skip_build:
        build_mcfost(args.compiler, args.jobs, dry_run=False)  # building is safe/reversible, always do it
    else:
        print("\n== Skipping build ==")

    if not args.skip_generate:
        generate_reference_data(dry_run=False)  # regenerating local test_data/ is safe/reversible too
    else:
        print("\n== Skipping reference data generation ==")

    if args.tarball:
        tarball = args.tarball.resolve()
        if not tarball.is_file():
            sys.exit(f"error: {tarball} does not exist")
        print(f"\n== Using existing tarball {tarball.name} ==")
    else:
        sha = git_sha()
        tag_date = date.today().strftime("%d_%m_%Y")
        # Packaging can produce a very large (hundreds of MB) file, so - unlike
        # build/generate above - it respects --dry-run too.
        tarball = make_tarball(sha, tag_date, dry_run=args.dry_run)

    if args.dry_run:
        if args.remote and not args.no_upload:
            upload_rclone(tarball, args.remote, dry_run=True)
            update_get_test_data_sh(tarball.name, "<new-file-id>", dry_run=True)
        print(
            "\nThis was a dry run: build/generate ran for real (if not skipped), "
            "but no tarball was written and nothing was uploaded or rewritten. "
            "Re-run without --dry-run to actually package and publish."
        )
        return

    if args.no_upload:
        print(f"\nDone. Tarball ready at {tarball} - not uploading (--no-upload).")
        return

    if not args.remote:
        print(
            "\nNo --remote given (and MCFOST_TEST_DATA_REMOTE is unset) - "
            "skipping upload. Tarball is ready at:", tarball
        )
        return

    # We only reach here without --dry-run, so this genuinely uploads and
    # rewrites get_test_data.sh.
    file_id, link = upload_rclone(tarball, args.remote, dry_run=False)
    if file_id:
        print(f"  uploaded, file id = {file_id}")
        print(f"  link: {link}")
        update_get_test_data_sh(tarball.name, file_id, dry_run=False)
        print("\nDone. Review the change to get_test_data.sh and commit it.")
    else:
        print(
            "\nUpload step did not produce a usable file ID - "
            "get_test_data.sh was left untouched. Update it by hand."
        )


if __name__ == "__main__":
    main()
