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

	ga := fem.NewGeneralAssembler(nodes, fem.DofPosX|fem.DofPosY)
	err := ga.AddIsoparametric(elements.Quad4{}, material.Axisymmetric(), 1, func(i int) (elem []int, xC, yC r3.Vec) {
		return elems[i][:], r3.Vec{}, r3.Vec{}
	})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("K=\n%.3g\n", lap.Formatted(ga.Ksolid()))
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
