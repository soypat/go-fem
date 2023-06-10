package fem_test

import (
	"fmt"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/elements"
	"github.com/soypat/lap"
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
	ga := fem.NewGeneralAssembler(nod, fem.DofPos)
	err := ga.AddIsoparametric3(elemType, steel, len(elems), func(i int) ([]int, r3.Vec, r3.Vec) {
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

// Calculates the stiffness matrix of a single isoparametric tetrahedron element.
func ExampleAnisotripicConstitutive() {
	carbonFiber := fem.TransverselyIsotropicMaterial{
		Ex:        235e9,
		Exy:       15e9,
		Gxy:       28e9,
		PoissonXY: 0.2,
		PoissonYZ: 0.25,
	}
	C, err := carbonFiber.Constitutive()
	if err != nil {
		panic(err)
	}
	fmt.Printf("C=\n%.3g", lap.Formatted(C))
	//Output:
	// C=
	// ⎡2.37e+11  4.03e+09  4.03e+09         0         0         0⎤
	// ⎢4.03e+09  1.61e+10  4.07e+09         0         0         0⎥
	// ⎢4.03e+09  4.07e+09  1.61e+10         0         0         0⎥
	// ⎢       0         0         0     6e+09         0         0⎥
	// ⎢       0         0         0         0   2.8e+10         0⎥
	// ⎣       0         0         0         0         0   2.8e+10⎦
}

func ExampleAxisymmetricConstitutive() {
	steel := fem.IsotropicMaterial{E: 1000, Poisson: 0.33}
	steelAxisymmetric := fem.AxisymmetricIsotropicMaterial{
		IsotropicMaterial: steel,
	}
	C, err := steelAxisymmetric.Constitutive()
	if err != nil {
		panic(err)
	}
	fmt.Printf("C=\n%.3g", lap.Formatted(C))
	//Output:

}
