package elements

import (
	"fmt"
	"math/rand"
	"testing"

	"github.com/soypat/go-fem"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/floats/scalar"
	"gonum.org/v1/gonum/spatial/r3"
)

var elements3 = []iso3{
	Tetra4{},
	Tetra10{},
	Hexa8{},
	Hexa20{},
}

type iso3 interface {
	fem.Isoparametric
	fmt.Stringer
	volume() float64
}

func TestBasis3d(t *testing.T) {
	const tol = 1e-16
	for _, element := range elements3 {
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

func TestBasisDiff3d(t *testing.T) {
	const (
		tol = 1e-11
		h   = 1e-4
	)
	ix := r3.Vec{X: h / 2}
	iy := r3.Vec{Y: h / 2}
	iz := r3.Vec{Z: h / 2}
	rng := rand.New(rand.NewSource(1))
	for _, element := range elements3 {
		t.Run(element.String(), func(t *testing.T) {
			for i := 0; i < 100; i++ {
				p := r3.Vec{X: rng.Float64(), Y: rng.Float64(), Z: rng.Float64()}
				// X calculation
				pxp := r3.Add(p, ix)
				pxm := r3.Sub(p, ix)
				ffxm := element.Basis(pxm)
				ffxp := element.Basis(pxp)
				// Differentiate and store in ffxp.
				floats.Sub(ffxp, ffxm)
				floats.Scale(1/h, ffxp)
				got := element.BasisDiff(p)
				if !floats.EqualApprox(got[:element.LenNodes()], ffxp, tol) {
					t.Error("X basis differential not equal", got, ffxp)
					return
				}

				// Y Differentiate
				pyp := r3.Add(p, iy)
				pym := r3.Sub(p, iy)
				ffym := element.Basis(pym)
				ffyp := element.Basis(pyp)
				floats.Sub(ffyp, ffym)
				floats.Scale(1/h, ffyp)
				if !floats.EqualApprox(got[element.LenNodes():element.LenNodes()*2], ffyp, tol) {
					t.Error("Y basis differential not equal", got, ffyp)
					return
				}

				// Z Differentiate
				pzp := r3.Add(p, iz)
				pzm := r3.Sub(p, iz)
				ffzm := element.Basis(pzm)
				ffzp := element.Basis(pzp)
				floats.Sub(ffzp, ffzm)
				floats.Scale(1/h, ffzp)
				if !floats.EqualApprox(got[2*element.LenNodes():], ffzp, tol) {
					t.Error("Y basis differential not equal", got, ffzp)
					return
				}
			}
		})
	}
}

func TestQuadrature3d(t *testing.T) {
	const (
		tol = 1e-11
	)
	for _, element := range elements3 {
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
			elementVolume := element.volume()
			if !scalar.EqualWithinAbs(sumW, elementVolume, tol) {
				t.Errorf("Sum of weights must be element volume %g, got %g", elementVolume, sumW)
			}
			if !scalar.EqualWithinAbs(sumN, elementVolume, tol) {
				t.Errorf("Sum of weights times shape functions must be element volume %g, got %g", elementVolume, sumW)
			}
		})
	}
}
