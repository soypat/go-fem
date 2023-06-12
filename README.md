# go-fem
Go packages for finite element analysis/methods.

#### Roadmap
* ~~3D Isoparametric element assembly.~~
* ~~Better sparse assembler~~ As fast as high performance libraries.
* ~~2D Isoparametric element assembly.~~
* Shell and plate element assembly.
* Arbitrary element assembly.
* Better define Element3 API.
* Add constraints assembly.
  - Lagrange Multipliers or Penalty method.
* Stiffness matrix utility functions for better conditioning and troubleshooting.
* Stress extraction from displacements.


#### Example
```go
package main

import (
	"fmt"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/elements"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

// Calculates the stiffness matrix of a single isoparametric tetrahedron element.
// and prints the result.
func main() {
	elemType := elements.Tetra4{}
	steel := fem.IsotropicMaterial{E: 200e9, Poisson: 0.3}
	nod := []r3.Vec{
		{X: 0, Y: 0, Z: 0},
		{X: 1}, // is (1,0,0) point
		{Y: 1},
		{Z: 1},
	}
	elems := [][]int{
		{0, 1, 2, 3}, // one single element.
	}
	ga := fem.NewGeneralAssembler(nod, fem.DofPos)
	err := ga.AddIsoparametric3(elemType, steel, len(elems), func(i int) ([]int, r3.Vec, r3.Vec) {
		return elems[i], r3.Vec{}, r3.Vec{}
	})
	if err != nil {
		panic(err)
	}
	fmt.Printf("K=\n%.5g", mat.Formatted(ga.Ksolid()))
}
```

**Output:**
```

K=
⎡ 4.2308e+11   1.9231e+11   1.9231e+11  -2.6923e+11  -7.6923e+10  -7.6923e+10  -7.6923e+10  -1.1538e+11            0  -7.6923e+10            0  -1.1538e+11⎤
⎢ 1.9231e+11   4.2308e+11   1.9231e+11  -1.1538e+11  -7.6923e+10            0  -7.6923e+10  -2.6923e+11  -7.6923e+10            0  -7.6923e+10  -1.1538e+11⎥
⎢ 1.9231e+11   1.9231e+11   4.2308e+11  -1.1538e+11            0  -7.6923e+10            0  -1.1538e+11  -7.6923e+10  -7.6923e+10  -7.6923e+10  -2.6923e+11⎥
⎢-2.6923e+11  -1.1538e+11  -1.1538e+11   2.6923e+11            0            0            0   1.1538e+11            0            0            0   1.1538e+11⎥
⎢-7.6923e+10  -7.6923e+10            0            0   7.6923e+10            0   7.6923e+10            0            0            0            0            0⎥
⎢-7.6923e+10            0  -7.6923e+10            0            0   7.6923e+10            0            0            0   7.6923e+10            0            0⎥
⎢-7.6923e+10  -7.6923e+10            0            0   7.6923e+10            0   7.6923e+10            0            0            0            0            0⎥
⎢-1.1538e+11  -2.6923e+11  -1.1538e+11   1.1538e+11            0            0            0   2.6923e+11            0            0            0   1.1538e+11⎥
⎢          0  -7.6923e+10  -7.6923e+10            0            0            0            0            0   7.6923e+10            0   7.6923e+10            0⎥
⎢-7.6923e+10            0  -7.6923e+10            0            0   7.6923e+10            0            0            0   7.6923e+10            0            0⎥
⎢          0  -7.6923e+10  -7.6923e+10            0            0            0            0            0   7.6923e+10            0   7.6923e+10            0⎥
⎣-1.1538e+11  -1.1538e+11  -2.6923e+11   1.1538e+11            0            0            0   1.1538e+11            0            0            0   2.6923e+11⎦
```

#### Benchmarks
Takes ~33s to assemble 680k dofs for 1.3 million 4 node 3D tetrahedrons on a low end AMD Zen2 CPU:
```
$ go test -bench=. -benchmem .
goos: linux
goarch: amd64
pkg: github.com/soypat/go-fem
cpu: AMD Ryzen 5 3400G with Radeon Vega Graphics    
BenchmarkTetra4Assembly/105_dofs,_48_elems-8                1218            932708 ns/op          155646 B/op        744 allocs/op
BenchmarkTetra4Assembly/3723_dofs,_5376_elems-8               10         109067373 ns/op        10928022 B/op      74500 allocs/op
BenchmarkTetra4Assembly/27027_dofs,_46080_elems-8              1        1023632489 ns/op        88685256 B/op     634443 allocs/op
BenchmarkTetra4Assembly/206115_dofs,_380928_elems-8            1        8849933029 ns/op        717490648 B/op   5229747 allocs/op
BenchmarkTetra4Assembly/684723_dofs,_1299456_elems-8           1        32997698741 ns/op       2742313160 B/op 17888792 allocs/op
```
