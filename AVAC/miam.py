from scipy.io import FortranFile


def readline(file):
    s = file.readline()
    return s[:s.rfind(" ")]

with open("_output/fgout0001.q0001") as file:
    grid_number = int(readline(file))
    AMR_level = int(readline(file))
    mx = int(readline(file))
    my = int(readline(file))
    xlow = float(readline(file))
    ylow = float(readline(file))
    dx = float(readline(file))
    dy = float(readline(file))

f = FortranFile( '_output/fgout0001.b0001', 'r' )
x = f.read_record( dtype=f'({my},{my})float32' )
print(x)
