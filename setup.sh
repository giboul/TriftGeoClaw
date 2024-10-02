sudo apt install libgdal-dev
pip install -r pyrequirements.txt
pip install --src=$CLAW/.. --no-build-isolation -e git+https://github.com/clawpack/clawpack.git@v5.11.0#egg=clawpack
echo "PETSC_DIR/PETSC_ARCH/lib/... could match:"
find / -name petsc 2>/dev/null
