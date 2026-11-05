#!/usr/bin/env bash
set -euxo pipefail

cargo install \
    --locked \
    --root "${PREFIX}" \
    --path . \
    --bin kun_peng

cargo-bundle-licenses \
    --format yaml \
    --output THIRDPARTY.yml
