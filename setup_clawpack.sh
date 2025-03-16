python3 -m venv claw
export CLAW=$PWD/claw/src/clawpack
echo "export CLAW=$CLAW" >> $PWD/claw/bin/activate
source claw/bin/activate

pip install meson-python
pip install ninja
pip install -r requirements.txt
pip install --no-build-isolation -e git+https://github.com/clawpack/clawpack.git@v5.11.0#egg=clawpack
sudo apt install liblapack-dev libopenblas-dev
