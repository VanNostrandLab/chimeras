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

cp "${cwd}"/source/*.py "${venv}"/bin/
sed -i "s|ENVIRONMENT|${venv}/environment.sh|" "${venv}/bin/chimeras.py"
sed -i "s|APP|${cwd}/chimeras|" "${venv}/bin/chimeras.py"

echo "# Successfully install and update chimeras"
echo "#"
echo "#    Run ${cwd}/chimeras -h to see the usage."
echo "#"
