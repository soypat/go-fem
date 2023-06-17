package fem_test

import (
	"fmt"
	"log"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/constitution/solids"
	"github.com/soypat/go-fem/elements"
	"github.com/soypat/lap"
	"gonum.org/v1/gonum/spatial/r3"
)

func ExampleGeneralAssembler_AddIsoparametric_axisymmetric() {
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
	//Output:
	// K=
	// ⎡ 2.93e+04   5.23e+03   1.45e+04   2.06e+03  -1.63e+04  -6.45e+03  -2.77e+04       -848⎤
	// ⎢ 5.23e+03   1.11e+05  -2.36e+03   6.14e+04  -7.37e+03  -6.19e+04        848  -1.11e+05⎥
	// ⎢ 1.45e+04  -2.36e+03    3.6e+04  -8.59e+03  -3.37e+04   3.58e+03  -1.63e+04   7.37e+03⎥
	// ⎢ 2.06e+03   6.14e+04  -8.59e+03   1.36e+05  -3.58e+03  -1.36e+05   6.45e+03  -6.19e+04⎥
	// ⎢-1.63e+04  -7.37e+03  -3.37e+04  -3.58e+03    3.6e+04   8.59e+03   1.45e+04   2.36e+03⎥
	// ⎢-6.45e+03  -6.19e+04   3.58e+03  -1.36e+05   8.59e+03   1.36e+05  -2.06e+03   6.14e+04⎥
	// ⎢-2.77e+04        848  -1.63e+04   6.45e+03   1.45e+04  -2.06e+03   2.93e+04  -5.23e+03⎥
	// ⎣     -848  -1.11e+05   7.37e+03  -6.19e+04   2.36e+03   6.14e+04  -5.23e+03   1.11e+05⎦
}

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

func ExampleGeneralAssembler_AddIsoparametric_axisymmetricThermal() {
	// We assemble a single 2D axisymmetric element.
	materialTherm := solids.IsotropicConductivity{K: 45} // 45 W/mK for steel.
	nodes := []r3.Vec{
		{X: 20, Y: 0},
		{X: 30, Y: 0},
		{X: 30, Y: 1},
		{X: 20, Y: 1},
	}
	elems := [][4]int{
		{0, 1, 2, 3},
	}
	elemType := elements.Quad4{NodeDofs: fem.DofPosX}
	ga := fem.NewGeneralAssembler(nodes, elemType.Dofs())
	err := ga.AddIsoparametric(elemType, materialTherm.Axisymmetric(), 1, func(i int) (elem []int, xC, yC r3.Vec) {
		return elems[i][:], r3.Vec{}, r3.Vec{}
	})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("K=\n%.3g\n", lap.Formatted(ga.Ksolid()))
	//Output:
	// K=
	// ⎡  3.4e+03   1.84e+03  -1.89e+03  -3.36e+03⎤
	// ⎢ 1.84e+03   4.18e+03   -4.1e+03  -1.89e+03⎥
	// ⎢-1.89e+03   -4.1e+03   4.18e+03   1.84e+03⎥
	// ⎣-3.36e+03  -1.89e+03   1.84e+03    3.4e+03⎦
}
