package fem_test

import (
	"fmt"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/constitution/solids"
	"github.com/soypat/go-fem/elements"
	"github.com/soypat/lap"
	"gonum.org/v1/gonum/spatial/r3"
)

// Calculates the stiffness matrix of a single isoparametric tetrahedron element.
func ExampleGeneralAssembler() {
	elemType := elements.Tetra4{}
	steel := solids.Isotropic{E: 200e9, Poisson: 0.3}
	nod := []r3.Vec{
		{X: 0, Y: 0, Z: 0},
		{X: 1}, // is (1,0,0) point
		{Y: 1},
		{Z: 1},
	}
	elems := [][]int{
		{0, 1, 2, 3}, // one single element.
	}
	ga := fem.NewGeneralAssembler(nod, elemType.Dofs())
	err := ga.AddIsoparametric(elemType, steel.Solid3D(), len(elems), func(i int) ([]int, r3.Vec, r3.Vec) {
		return elems[i], r3.Vec{}, r3.Vec{}
	})
	if err != nil {
		panic(err)
	}
	fmt.Printf("K=\n%.5g", lap.Formatted(ga.Ksolid()))
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

func ExampleGeneralAssembler_AddIsoparametric_quad4PlaneStress() {
	// We assemble a single Quad4 element with a unitary isotropic material.
	material := solids.Isotropic{E: 1, Poisson: 0.3}
	nodes := []r3.Vec{
		{X: 0, Y: 0},
		{X: 1, Y: 0},
		{X: 1, Y: 1},
		{X: 0, Y: 1},
	}
	elemType := elements.Quad4{}
	elems := [][]int{
		{0, 1, 2, 3},
	}

	ga := fem.NewGeneralAssembler(nodes, elemType.Dofs())
	err := ga.AddIsoparametric(elemType, material.PlaneStess(), len(elems), func(i int) (elem []int, xC, yC r3.Vec) {
		return elems[i], r3.Vec{}, r3.Vec{}
	})
	if err != nil {
		panic(err)
	}
	fmt.Printf("K=\n%.3f", lap.Formatted(ga.Ksolid()))
	//Output:
	// K=
	// ⎡ 0.495   0.179  -0.302  -0.014  -0.247  -0.179   0.055   0.014⎤
	// ⎢ 0.179   0.495   0.014   0.055  -0.179  -0.247  -0.014  -0.302⎥
	// ⎢-0.302   0.014   0.495  -0.179   0.055  -0.014  -0.247   0.179⎥
	// ⎢-0.014   0.055  -0.179   0.495   0.014  -0.302   0.179  -0.247⎥
	// ⎢-0.247  -0.179   0.055   0.014   0.495   0.179  -0.302  -0.014⎥
	// ⎢-0.179  -0.247  -0.014  -0.302   0.179   0.495   0.014   0.055⎥
	// ⎢ 0.055  -0.014  -0.247   0.179  -0.302   0.014   0.495  -0.179⎥
	// ⎣ 0.014  -0.302   0.179  -0.247  -0.014   0.055  -0.179   0.495⎦
}
