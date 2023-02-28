package elements

import (
	"errors"

	"github.com/soypat/go-fem"
	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

// Beam2dof6 represents a 1D beam element with 6 dof on each of it's two nodes.
// The beam is by default oriented in X direction.
type Beam2dof6 struct {
	// Cross-sectional area of beam (YZ plane in local coordinates)
	A float64
	// Area moment of inertia of beam corresponding to couple on local Y axis.
	Iy float64
	// Area moment of inertia of beam corresponding to couple on local Z axis.
	Iz float64
	// Longitudinal moment of inertia of beam corresponding to torsion on local X axis.
	J float64
	e float64
	g float64
}

var _ fem.Element3 = (*Beam2dof6)(nil)

func (b *Beam2dof6) LenNodes() int { return 2 }

func (b *Beam2dof6) Dofs() fem.DofsFlag { return fem.Dof6 }

func (b *Beam2dof6) CopyK(dst *mat.Dense, v []r3.Vec) error {
	if len(v) != 2 {
		return errors.New("need 2 nodes")
	}
	L := r3.Norm(r3.Sub(v[0], v[1]))
	E := b.e
	G := b.g
	X := b.A * E / L
	Y4 := 2 * E * b.Iz / L
	Y3 := Y4 * 2
	Y2 := 3 * Y4 / L
	Y1 := 2 * Y2 / L
	Z4 := 2 * E * b.Iy / L
	Z3 := Z4 * 2
	Z2 := 3 * Z4 / L
	Z1 := 2 * Z2 / L
	S := G * b.J / L
	dst.SetRawMatrix(blas64.General{
		Rows:   12,
		Cols:   12,
		Stride: 12,
		Data: []float64{
			X, 0, 0, 0, 0, 0, -X, 0, 0, 0, 0, 0,
			0, Y1, 0, 0, 0, Y2, 0, -Y1, 0, 0, 0, Y2,
			0, 0, Z1, 0, -Z2, 0, 0, 0, -Z1, 0, -Z2, 0,
			0, 0, 0, S, 0, 0, 0, 0, 0, -S, 0, 0,
			0, 0, -Z2, 0, Z3, 0, 0, 0, Z2, 0, Z4, 0,
			0, Y2, 0, 0, 0, Y3, 0, -Y2, 0, 0, 0, Y4,
			-X, 0, 0, 0, 0, 0, X, 0, 0, 0, 0, 0,
			0, -Y1, 0, 0, 0, -Y2, 0, Y1, 0, 0, 0, -Y2,
			0, 0, -Z1, 0, Z2, 0, 0, 0, Z1, 0, Z2, 0,
			0, 0, 0, -S, 0, 0, 0, 0, 0, S, 0, 0,
			0, 0, -Z2, 0, Z4, 0, 0, 0, Z2, 0, Z3, 0,
			0, Y2, 0, 0, 0, Y4, 0, -Y2, 0, 0, 0, Y3,
		},
	})
	return nil
}

func (b *Beam2dof6) SetConstitutive(c fem.Constituter) error {
	sp, err := extractSolidProps(c)
	if err != nil {
		return errors.New("invalid contitutive parameters: " + err.Error())
	}
	// TODO(soypat): Better way to construct beam than with standard two values?
	// Probably incorrect for non-isotropic materials
	b.e = (sp.Ex + sp.Ey + sp.Ez) / 3
	b.g = (sp.Gxy + sp.Gxz + sp.Gyz) / 3
	return nil
}

type solidProperties struct {
	// Young modulii.
	Ex, Ey, Ez float64
	// Supposes Vxy == Vyx etc (symmetric constitutive matrix).
	Vyz, Vxz, Vxy float64
	// Shear modulii.
	Gyz, Gxz, Gxy float64
}

func extractSolidProps(c fem.Constituter) (solidProperties, error) {
	var sp solidProperties
	s := mat.NewDense(6, 6, nil)
	C, err := c.Constitutive()
	if err != nil {
		return solidProperties{}, err
	}
	err = s.Inverse(C)
	if err != nil {
		return solidProperties{}, err
	}
	// Extracting data from constitutive matrix.
	// http://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_no_solutions.pdf
	sp.Ex = 1 / s.At(0, 0)
	sp.Ey = 1 / s.At(1, 1)
	sp.Ez = 1 / s.At(2, 2)
	sp.Vxy = -sp.Ey * s.At(0, 1)
	sp.Vxz = -sp.Ez * s.At(0, 2)
	sp.Vyz = -sp.Ez * s.At(1, 2)
	sp.Gyz = 1 / s.At(3, 3)
	sp.Gxz = 1 / s.At(4, 4)
	sp.Gxy = 1 / s.At(5, 5)
	return sp, nil
}
