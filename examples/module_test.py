"""
Version: 
    1.0

Copyright Notice:
    Copyright (c) 2021, Claudia Campos
    All rights reserved.

License Notice:

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

from parallelcoordinates import ParallelCoordinates
from pandas import read_csv
import numpy as np

data = np.random.laplace(loc=1, scale=5.0, size=(100,7))

g = ParallelCoordinates(data)

input("Press Enter to continue...")

g.updateAlpha(0.25)

input("Press Enter to continue...")

g.updateColor('#5900ff')

input("Press Enter to continue...")

g.setScale(1.2)

g.setAutoSpace(False)

g.setPosition(np.random.uniform(-0.5, g.nDims - 0.5, g.nDims))

input("Press Enter to continue...")

g.setPermutation(np.random.permutation(g.nDims))

input("Press Enter to continue...")

g.resetAttributes()

print("Axis can be dragged using the left mouse key.")
print("Zoom in and zoom out using the mouse scroll.")