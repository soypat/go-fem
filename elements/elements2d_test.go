package elements

import (
	"fmt"
	"math/rand"
	"testing"

	"github.com/soypat/go-fem"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/floats/scalar"
	"gonum.org/v1/gonum/spatial/r2"
)

var elements2 = []iso2{
	Quad8{},
	Quad4{},
	Triangle3{},
	// Triangle6{},
}

type iso2 interface {
	fem.Isoparametric2
	fmt.Stringer
	area() float64
}

func TestBasis2(t *testing.T) {
	const tol = 1e-16
	for _, element := range elements2 {
		t.Run(element.String(), func(t *testing.T) {
			nod := element.IsoparametricNodes()
			for i, n := range nod {
				ff := element.Basis(n)
				for j := range ff {
					if i == j && !scalar.EqualWithinAbs(ff[j], 1, tol) {
						t.Errorf("basis%d must be 1 at corresponding node, got %f", i, ff[j])
					} else if i != j && !scalar.EqualWithinAbs(ff[j], 0, tol) {
						t.Errorf("basis%d must be 0 at non corresponding nodes, got %f", i, ff[j])
					}
				}
			}
		})
	}
}

func TestBasisDiff2d(t *testing.T) {
	const (
		tol = 1e-11
		h   = 1e-4
	)
	ix := r2.Vec{X: h / 2}
	iy := r2.Vec{Y: h / 2}
	rng := rand.New(rand.NewSource(1))
	for _, element := range elements2 {
		t.Run(element.String(), func(t *testing.T) {
			for i := 0; i < 100; i++ {
				p := r2.Vec{X: rng.Float64(), Y: rng.Float64()}
				// X calculation
				pxp := r2.Add(p, ix)
				pxm := r2.Sub(p, ix)
				ffxm := element.Basis(pxm)
				ffxp := element.Basis(pxp)
				// Differentiate and store in ffxp.
				floats.Sub(ffxp, ffxm)
				floats.Scale(1/h, ffxp)
				got := element.BasisDiff(p)
				if !floats.EqualApprox(got[:element.LenNodes()], ffxp, tol) {
					t.Errorf("X basis differential not equal\ngot:%v\nwant:%v", got, ffxp)
					return
				}

				// Y Differentiate
				pyp := r2.Add(p, iy)
				pym := r2.Sub(p, iy)
				ffym := element.Basis(pym)
				ffyp := element.Basis(pyp)
				floats.Sub(ffyp, ffym)
				floats.Scale(1/h, ffyp)
				if !floats.EqualApprox(got[element.LenNodes():element.LenNodes()*2], ffyp, tol) {
					t.Errorf("Y basis differential not equal\ngot:%v\nwant:%v", got, ffyp)
					return
				}
			}
		})
	}
}

func TestQuadrature2d(t *testing.T) {
	const (
		tol = 1e-11
	)
	elements2WithQuads := []iso2{
		Quad4{QuadratureOrder: 1},
		Quad8{QuadratureOrder: 2},
		Quad8{QuadratureOrder: 1},
		Triangle3{QuadratureOrder: 2},
		Triangle6{QuadratureOrder: 1},
		Triangle6{QuadratureOrder: 3},
	}
	for _, element := range append(elements2, elements2WithQuads...) {
		t.Run(element.String(), func(t *testing.T) {
			wPositions, weights := element.Quadrature()
			if len(wPositions) != len(weights) {
				t.Fatal("Quadrature points and weights must have same length")
			}
			sumW := 0.0
			sumN := 0.0
			for i, w := range weights {
				sumW += w
				pos := wPositions[i]
				Nwpg := element.Basis(pos)
				for _, N := range Nwpg {
					sumN += N * w
				}
			}

			elementArea := element.area()
			if !scalar.EqualWithinAbs(sumW, elementArea, tol) {
				t.Errorf("Sum of weights must be element area %g, got %g", elementArea, sumW)
			}
			if !scalar.EqualWithinAbs(sumN, elementArea, tol) {
				t.Errorf("Sum of weights times shape functions must be element area %g, got %g", elementArea, sumW)
			}
		})
	}
}
