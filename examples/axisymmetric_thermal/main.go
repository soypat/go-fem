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
}
