#!/usr/bin/env bash
set -e

cwd="$( cd .; pwd -P )"
venv="${cwd}/venv"

read -r -d '' environment << EOF || true
#!/usr/bin/env bash

unset PYTHONPATH
export PATH="${venv}/bin:\$PATH"
EOF
echo "${environment}" > "${venv}/environment.sh"

read -r -d '' script << EOF || true
#!/usr/bin/env bash

source ${venv}/environment.sh
python ${venv}/bin/chimeras.py \$@
EOF
echo "${script}" > "${cwd}/chimeras"
chmod +x "${cwd}/chimeras"

source_py="${cwd}/source/chimeras.py"
bin_py="${venv}/bin/chimeras.py"
sed "s|ENVIRONMENT|${venv}/environment.sh|" "${source_py}" > "${bin_py}"
echo "# Successfully install and update chimeras"
echo "#"
echo "#    Run ${cwd}/chimeras -h to see the usage."
echo "#"
