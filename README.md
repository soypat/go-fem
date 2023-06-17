# go-fem
Go packages for finite element analysis/methods.

#### Roadmap
* ~~3D Isoparametric element assembly.~~
* ~~Better sparse assembler~~ As fast as high performance libraries.
* ~~2D Isoparametric element assembly.~~
* ~~Arbitrary Isoparametric element assembly.~~
* Stress extraction from displacements.
* Shell and plate element assembly.
* Better define Element3 API.
* Add constraints assembly.
  - Lagrange Multipliers or Penalty method.
* Stiffness matrix utility functions for better conditioning and troubleshooting.



#### Axisymmetric assembly example
```go
package main

import (
	"fmt"
	"log"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/constitution/solids"
	"github.com/soypat/go-fem/elements"
	"github.com/soypat/lap"
	"gonum.org/v1/gonum/spatial/r3"
)

func main() {
	// We assemble a single 2D axisymmetric element.
	// With a young modulus of 1000 and a poisson ratio of 0.33.
	// The units are undefined. If using SI, E would be 1kPa and node positions
	// would be in meters.
	material := solids.Isotropic{E: 1000, Poisson: 0.33}
	nodes := []r3.Vec{
		{X: 20, Y: 0},
		{X: 30, Y: 0},
		{X: 30, Y: 1},
		{X: 20, Y: 1},
	}
	elems := [][4]int{
		{0, 1, 2, 3},
	}
	elemType := elements.Quad4{}

	ga := fem.NewGeneralAssembler(nodes, elemType.Dofs())
	err := ga.AddIsoparametric(elemType, material.Axisymmetric(), 1, func(i int) (elem []int, xC, yC r3.Vec) {
		return elems[i][:], r3.Vec{}, r3.Vec{}
	})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("K=\n%.3g\n", lap.Formatted(ga.Ksolid()))
}
```

**Output:**
```
K=
⎡ 2.93e+04   5.23e+03   1.45e+04   2.06e+03  -1.63e+04  -6.45e+03  -2.77e+04       -848⎤
⎢ 5.23e+03   1.11e+05  -2.36e+03   6.14e+04  -7.37e+03  -6.19e+04        848  -1.11e+05⎥
⎢ 1.45e+04  -2.36e+03    3.6e+04  -8.59e+03  -3.37e+04   3.58e+03  -1.63e+04   7.37e+03⎥
⎢ 2.06e+03   6.14e+04  -8.59e+03   1.36e+05  -3.58e+03  -1.36e+05   6.45e+03  -6.19e+04⎥
⎢-1.63e+04  -7.37e+03  -3.37e+04  -3.58e+03    3.6e+04   8.59e+03   1.45e+04   2.36e+03⎥
⎢-6.45e+03  -6.19e+04   3.58e+03  -1.36e+05   8.59e+03   1.36e+05  -2.06e+03   6.14e+04⎥
⎢-2.77e+04        848  -1.63e+04   6.45e+03   1.45e+04  -2.06e+03   2.93e+04  -5.23e+03⎥
⎣     -848  -1.11e+05   7.37e+03  -6.19e+04   2.36e+03   6.14e+04  -5.23e+03   1.11e+05⎦
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
