#!/usr/bin/env bash
set -e -o pipefail
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SRC_DIR="${SCRIPT_DIR}/IsisRNA_ToolsSrc"
cd "${SRC_DIR}"
./configure --prefix="${SCRIPT_DIR}/IsisRNA_Tools"
make -j 4 release
cd  "${SCRIPT_DIR}"
echo "Deleting ${SRC_DIR}!!"
rm -r -f "${SRC_DIR}"
rm -r -f ${BASH_SOURCE[0]}
