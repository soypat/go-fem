package main

import (
	"fmt"
	"log"

	"github.com/soypat/go-fem"
	"github.com/soypat/go-fem/constitution/solids"
	"github.com/soypat/go-fem/elements"
	"github.com/soypat/lap"
	"gonum.org/v1/gonum/mat"
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
	// var jac, dNxy mat.Dense
	// upg, wpg := elemType.Quadrature()
	// err = ga.IntegrateIsoparametric(elemType, 1, func(i int) (elem []int, xC, yC r3.Vec) {
	// 	return elems[i][:], r3.Vec{}, r3.Vec{}
	// }, func(elemIdx int, elemNodes []float64, elemDofs []int) error {
	// 	forceElem := lap.NewDenseVector(len(elemDofs), nil)

	// 	for ipg := range dNpg {
	// 		N := Npg[ipg]
	// 		jac.Mul(dNpg[ipg], elemNodes)
	// 		dNxy.Solve(&jac, dNpg[ipg])
	// 		radius := lap.Dot(N, elemNodes.ColView(0))
	// 		forceElem.AddVec(forceElem, N)
	// 		fmt.Printf("J%d%d=\n%.3g\n", elemIdx, ipg, lap.Formatted(&jac))
	// 	}
	// 	return nil
	// })
	// var forces lap.Vector
	// forces = lap.NewDenseVector(ga.TotalDofs(), []float64{1.3750 * 1e3, 0 * 1e3, 1.7917 * 1e3, 0 * 1e3, 1.7917 * 1e3, 0 * 1e3, 1.3750 * 1e3, 0})
	// fixedDofs := []int{1}
	// forces = lap.SliceExcludeVec(forces, fixedDofs)

	// var K lap.Matrix
	// K = ga.Ksolid()
	// K = lap.SliceExclude(K, fixedDofs, fixedDofs)
	// var displacementsReduced mat.VecDense
	// err = displacementsReduced.SolveVec(lapmat{K}, lapvec{forces})
	// if err != nil {
	// 	log.Fatal(err)
	// }
	// displacements := lap.NewDenseVector(ga.TotalDofs(), nil)
	// copied := lap.SliceExcludeVec(displacements, fixedDofs).CopyVec(&displacementsReduced)
	// fmt.Println("displacements", copied, lap.Formatted(displacements))
}

type lapvec struct {
	lap.Vector
}

func (v lapvec) T() mat.Matrix {
	return lapmat{lap.T(v.Vector)}
}

type lapmat struct {
	lap.Matrix
}

func (m lapmat) T() mat.Matrix {
	if ter, ok := m.Matrix.(mat.Matrix); ok {
		return ter.T()
	}
	return lapmat{lap.T(m.Matrix)}
}

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
