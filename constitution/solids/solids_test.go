package solids_test

import (
	"fmt"

	"github.com/soypat/go-fem/constitution/solids"
	"github.com/soypat/lap"
)

// Calculates the stiffness matrix of a single isoparametric tetrahedron element.
func ExampleTransverselyIsotropic_carbonFiber() {
	carbonFiber := solids.TransverselyIsotropic{
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
