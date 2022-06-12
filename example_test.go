package fem_test

import (
	"fmt"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/elements"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

// Calculates the stiffness matrix of a single isoparametric tetrahedron element.
func ExampleGeneralAssembler() {
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
	ga := fem.NewGeneralAssembler(nod, fem.DofU)
	err := ga.AddIsoparametric3s(elemType, steel, len(elems), func(i int) []int {
		return elems[i]
	})
	if err != nil {
		panic(err)
	}
	fmt.Printf("K=\n%.5g", mat.Formatted(ga.Ksolid()))
	// Output:
	//  K=
	// ⎡ 4.2308e+11   1.9231e+11   1.9231e+11  -2.6923e+11  -7.6923e+10  -7.6923e+10  -7.6923e+10  -1.1538e+11            0  -7.6923e+10            0  -1.1538e+11⎤
	// ⎢ 1.9231e+11   4.2308e+11   1.9231e+11  -1.1538e+11  -7.6923e+10            0  -7.6923e+10  -2.6923e+11  -7.6923e+10            0  -7.6923e+10  -1.1538e+11⎥
	// ⎢ 1.9231e+11   1.9231e+11   4.2308e+11  -1.1538e+11            0  -7.6923e+10            0  -1.1538e+11  -7.6923e+10  -7.6923e+10  -7.6923e+10  -2.6923e+11⎥
	// ⎢-2.6923e+11  -1.1538e+11  -1.1538e+11   2.6923e+11            0            0            0   1.1538e+11            0            0            0   1.1538e+11⎥
	// ⎢-7.6923e+10  -7.6923e+10            0            0   7.6923e+10            0   7.6923e+10            0            0            0            0            0⎥
	// ⎢-7.6923e+10            0  -7.6923e+10            0            0   7.6923e+10            0            0            0   7.6923e+10            0            0⎥
	// ⎢-7.6923e+10  -7.6923e+10            0            0   7.6923e+10            0   7.6923e+10            0            0            0            0            0⎥
	// ⎢-1.1538e+11  -2.6923e+11  -1.1538e+11   1.1538e+11            0            0            0   2.6923e+11            0            0            0   1.1538e+11⎥
	// ⎢          0  -7.6923e+10  -7.6923e+10            0            0            0            0            0   7.6923e+10            0   7.6923e+10            0⎥
	// ⎢-7.6923e+10            0  -7.6923e+10            0            0   7.6923e+10            0            0            0   7.6923e+10            0            0⎥
	// ⎢          0  -7.6923e+10  -7.6923e+10            0            0            0            0            0   7.6923e+10            0   7.6923e+10            0⎥
	// ⎣-1.1538e+11  -1.1538e+11  -2.6923e+11   1.1538e+11            0            0            0   1.1538e+11            0            0            0   2.6923e+11⎦
}
