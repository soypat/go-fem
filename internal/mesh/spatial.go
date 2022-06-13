package mesh

import (
	"math"

	"gonum.org/v1/gonum/spatial/r3"
)

func boxSize(b r3.Box) r3.Vec {
	return r3.Sub(b.Max, b.Min)
}

func boxVertices(a r3.Box) []r3.Vec {
	v := make([]r3.Vec, 8)
	v[0] = a.Min
	v[1] = r3.Vec{X: a.Min.X, Y: a.Min.Y, Z: a.Max.Z}
	v[2] = r3.Vec{X: a.Min.X, Y: a.Max.Y, Z: a.Min.Z}
	v[3] = r3.Vec{X: a.Min.X, Y: a.Max.Y, Z: a.Max.Z}
	v[4] = r3.Vec{X: a.Max.X, Y: a.Min.Y, Z: a.Min.Z}
	v[5] = r3.Vec{X: a.Max.X, Y: a.Min.Y, Z: a.Max.Z}
	v[6] = r3.Vec{X: a.Max.X, Y: a.Max.Y, Z: a.Min.Z}
	v[7] = a.Max
	return v
}

func boxCenter(a r3.Box) r3.Vec {
	return r3.Add(a.Min, r3.Scale(0.5, boxSize(a)))
}

// CenteredBox creates a Box with a given center and size.
// Negative components of size will be interpreted as zero.
func centeredBox(center, size r3.Vec) r3.Box {
	size = maxElem(size, r3.Vec{}) // set negative values to zero.
	half := r3.Scale(0.5, size)
	return r3.Box{Min: r3.Sub(center, half), Max: r3.Add(center, half)}
}

// MaxElem return a vector with the maximum components of two vectors.
func maxElem(a, b r3.Vec) r3.Vec {
	return r3.Vec{X: math.Max(a.X, b.X), Y: math.Max(a.Y, b.Y), Z: math.Max(a.Z, b.Z)}
}

// gradient returns the gradient of a scalar field. This also returns the normal
// vector of a sdf surface.
func gradient(p r3.Vec, tol float64, f func(r3.Vec) float64) r3.Vec {
	return r3.Vec{
		X: f(r3.Add(p, r3.Vec{X: tol})) - f(r3.Add(p, r3.Vec{X: -tol})),
		Y: f(r3.Add(p, r3.Vec{Y: tol})) - f(r3.Add(p, r3.Vec{Y: -tol})),
		Z: f(r3.Add(p, r3.Vec{Z: tol})) - f(r3.Add(p, r3.Vec{Z: -tol})),
	}
}

func divergence(p r3.Vec, tol float64, f func(r3.Vec) r3.Vec) float64 {
	dx := r3.Sub(f(r3.Add(p, r3.Vec{X: tol})), f(r3.Add(p, r3.Vec{X: -tol})))
	dy := r3.Sub(f(r3.Add(p, r3.Vec{Y: tol})), f(r3.Add(p, r3.Vec{Y: -tol})))
	dz := r3.Sub(f(r3.Add(p, r3.Vec{Z: tol})), f(r3.Add(p, r3.Vec{Z: -tol})))
	return dx.X + dy.Y + dz.Z
}

func laplacian(p r3.Vec, tol float64, f func(r3.Vec) float64) float64 {
	return divergence(p, tol, func(v r3.Vec) r3.Vec { return gradient(p, tol, f) })
}
