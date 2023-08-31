from reservoirflow import grids
import time


t1 = time.time()
grid = grids.Cartesian(nx=3, ny=1, nz=1, dx=10, dy=10, dz=5, phi=0.4, kx=10)
for i in range(1000):
    id = grid.get_cell_id((1, 0, 0))
    id = grid.get_cell_id((2, 0, 0))

print(id)
t2 = time.time()
t = (t2 - t1) / 60
print(f"finish in {t} minutes")
