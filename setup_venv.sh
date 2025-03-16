read CLAW
python3 -m venv claw
echo "export CLAW=$CLAW" >> $PWD/claw/bin/activate
source claw/bin/activate

pip install meson-python
pip install ninja
pip install -r requirements.txt
sudo apt install liblapack-dev libopenblas-dev
