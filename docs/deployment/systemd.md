# Running as a systemd Service

For long-running deployments (lab workstation, dedicated server) we
ship a **systemd user unit** that wraps `omnipath-metabo serve` so the
API restarts on failure and starts automatically at boot. This is what
we use for the public service at `metabo.omnipathdb.org` and is also
the recommended setup for self-hosted instances &mdash; for example,
when serving a Parquet cache built from your own COSMOS PKN variant.

The unit file is part of the repository at
[`deploy/systemd/omnipath-metabo.service`][unit]. It assumes:

- the repository is checked out at `~/dev/omnipath-metabo/`,
- a uv-managed virtualenv exists at `~/dev/omnipath-metabo/.venv/`
  (created by `uv sync` or `pip install -e .[server]`),
- the Parquet cache (built once via the Python or `cosmos-pkn` CLI
  workflows in [Quickstart](../quickstart.md)) lives at the default
  location `~/.cache/omnipath-metabo/` or wherever
  `OMNIPATH_METABO_CACHE_DIR` points.

[unit]: https://github.com/saezlab/omnipath-metabo/blob/main/deploy/systemd/omnipath-metabo.service

## Install

Run the four commands below as the user that will own the service
(does **not** need to be root):

```bash
git clone https://github.com/saezlab/omnipath-metabo.git ~/dev/omnipath-metabo
mkdir -p ~/.config/systemd/user
ln -sf ~/dev/omnipath-metabo/deploy/systemd/omnipath-metabo.service \
       ~/.config/systemd/user/omnipath-metabo.service
systemctl --user daemon-reload
systemctl --user enable --now omnipath-metabo.service
```

To survive reboots without an active login session, enable lingering
once as root:

```bash
sudo loginctl enable-linger "$USER"
```

## Configure

The unit defines two environment variables with sensible defaults:

```ini
Environment=OMNIPATH_METABO_HOST=127.0.0.1
Environment=OMNIPATH_METABO_PORT=8084
```

Override them &mdash; or add `OMNIPATH_METABO_CACHE_DIR` to point at a
custom Parquet cache &mdash; via a drop-in instead of editing the
shipped file:

```bash
systemctl --user edit omnipath-metabo
```

```ini
[Service]
Environment=OMNIPATH_METABO_PORT=9000
Environment=OMNIPATH_METABO_CACHE_DIR=/srv/metabo-cache-custom
```

If your environment needs additional setup (e.g. activating a Conda
profile or exporting an `LD_LIBRARY_PATH` on NixOS for numpy/pandas
C-extensions), drop the export statements into `~/dev/.envrc`. The unit
sources that file when present and ignores it otherwise.

## Operate

```bash
systemctl --user status  omnipath-metabo
systemctl --user restart omnipath-metabo
systemctl --user stop    omnipath-metabo
journalctl --user -u omnipath-metabo -f       # live logs
```

The service binds to `127.0.0.1` by default. To expose it publicly,
front it with nginx, Caddy, or an SSH reverse tunnel and terminate TLS
there &mdash; do not bind `0.0.0.0` directly on an untrusted network.

## System-wide alternative

If you prefer a system service (e.g. a dedicated `omnipath` system
account with no login shell), copy the same unit to
`/etc/systemd/system/omnipath-metabo.service`, replace `%h` with the
target home directory, add `User=` / `Group=` lines under `[Service]`,
and run `systemctl daemon-reload && systemctl enable --now
omnipath-metabo.service` as root. No other changes are required.

## Multiple instances on one host

For shared development hosts where several instances of the service
need to run side-by-side (e.g. a `staging` deployment alongside a
production deployment, or a per-feature dev URL), a templated unit is
shipped at [`deploy/systemd/omnipath-metabo@.service`][template].

Lay each instance out at `~/instances/<name>/` &mdash; with `src/`,
`.venv/`, `.env`, and a `.cache/` directory (often a symlink to a
shared parquet cache) &mdash; then enable:

```bash
systemctl --user enable --now omnipath-metabo@staging.service
```

Each instance gets its own port and journal stream; cache and code can
be shared with prod or kept separate per instance. The same
`loginctl enable-linger` once-per-user setup applies.

[template]: https://github.com/saezlab/omnipath-metabo/blob/main/deploy/systemd/omnipath-metabo@.service
