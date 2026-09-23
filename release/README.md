# Publishing this candidate

The standalone source folder is ready for local testing and review. A fresh
Git repository on branch `main` is initialized, with no remote or commits.
The research repository and its history are separate.

## Finish metadata

1. Choose the repository owner/name and provide its GitHub URL.
2. Supply the software author list and choose a license. Add the approved
   license text as `LICENSE`; do not treat the citation template as a license.
3. Fill `release/CITATION.cff.template` and save it as root `CITATION.cff`.
4. Update `release/metadata.json`, README's candidate notice, `VERSION`, and
   `CHANGELOG.md`. A paper DOI can be added later.

The archive builder automatically includes root LICENSE and CITATION.cff
when those files have been supplied. No license has been selected on the
owner's behalf.

## Create the GitHub repository

Create an empty public repository. Leave GitHub's initial README, license,
and .gitignore options unchecked because these will come from this folder.
Send the resulting repository URL back so that it can be connected.

After metadata review, the normal first-push sequence from this folder is:

```bash
git add .
git diff --cached --stat
git commit -m "Prepare initial numerical solver release"
git remote add origin https://github.com/OWNER/REPOSITORY.git
git push -u origin main
```

These commands are instructions, not actions already performed. Replace the
URL with the actual repository. Git author name/email must be configured by
the owner before committing. Review GitHub Actions results after the push.

## Share a stable version

Once the candidate is accepted, set the chosen release version, create its Git
tag, and create a GitHub Release with the corresponding notes. Build a fresh
source archive from that reviewed version:

```bash
python scripts/build_solver_release.py --output dist/rmi-solver-v0.1.0.zip
```

The builder refuses to overwrite existing archives. Attach the ZIP and its
`.sha256` sidecar to the GitHub Release. Share the release URL and ask a group
member to follow README setup and run the tests on another machine.

The current workflow follows GitHub's Python testing guide:
https://docs.github.com/en/actions/tutorials/build-and-test-code/python
It runs tests and examples; it has not yet run on GitHub and does not publish.
