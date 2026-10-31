# Releasing

For maintainers. Releases go to PyPI through
[Trusted Publishing](https://docs.pypi.org/trusted-publishers/): the
[`release` workflow](https://github.com/elmasnuryilmaz/citeqc/blob/main/.github/workflows/release.yml)
authenticates with a short-lived OpenID Connect token, so no API token is
stored in the repository.

| Trigger | What happens |
|---|---|
| push a tag `vX.Y.Z` | build, `twine check`, check that the tag equals `citeqc.__version__`, smoke-test the wheel, publish to **PyPI** |
| Actions, release, *Run workflow* | the same build and checks, then publish to **TestPyPI** (a dry run) |

## One-time setup

### 1. Register the trusted publishers

The project does not exist on PyPI yet, so it is registered as a *pending*
publisher; PyPI creates the project on the first successful upload. A pending
publisher does not reserve the name, so publish soon after setting it up.

1. Sign in at [pypi.org](https://pypi.org) (two-factor authentication must be
   on) and open **Account settings, Publishing**
   ([pypi.org/manage/account/publishing](https://pypi.org/manage/account/publishing/)).
2. Under **Add a new pending publisher**, choose **GitHub** and enter:

    | Field | Value |
    |---|---|
    | PyPI Project Name | `citeqc` |
    | Owner | `elmasnuryilmaz` |
    | Repository name | `citeqc` |
    | Workflow name | `release.yml` |
    | Environment name | `pypi` |

3. For dry runs, repeat on [test.pypi.org](https://test.pypi.org/manage/account/publishing/)
   (a separate account) with environment name `testpypi`.

### 2. Create the GitHub environments

In the repository, **Settings, Environments, New environment**: create `pypi`
and, if you use TestPyPI, `testpypi`. Recommended for `pypi`:

- **Required reviewers**: yourself, so every upload waits for your approval;
- **Deployment branches and tags**: only tags matching `v*`.

### 3. Enable GitHub Pages

**Settings, Pages, Build and deployment, Source: GitHub Actions.** The
[`docs` workflow](https://github.com/elmasnuryilmaz/citeqc/blob/main/.github/workflows/docs.yml)
then publishes the site on every push to `main`, at
<https://elmasnuryilmaz.github.io/citeqc/>. Pull requests build the site with
`mkdocs build --strict` but do not deploy it.

## Each release

1. Update the version in `src/citeqc/__init__.py` (`__version__`, the only
   place it is set), move the *Unreleased* notes in `CHANGELOG.md` under the new
   version with today's date, and update `version` and `date-released` in
   `CITATION.cff`.
2. Merge to `main` and wait for the `tests` workflow to pass.
3. Optional dry run: **Actions, release, Run workflow** on `main`, then

    ```bash
    pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ citeqc
    ```

    TestPyPI never accepts the same version twice; re-runs skip files that
    already exist.

4. Tag and push:

    ```bash
    git tag -a v0.2.0 -m "citeqc 0.2.0"
    ```

    ```bash
    git push origin v0.2.0
    ```

5. Create a GitHub release from the tag, with the changelog entry as notes. If
   the repository is connected to [Zenodo](https://zenodo.org), this also
   archives the release and mints a DOI you can cite.
