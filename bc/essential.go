package bc

import "github.com/soypat/go-fem"

type Essential struct {
	Nodes     []fem.DofsFlag
	Augmented []bool
	D         fem.DofsFlag
	AreFree   bool
	dofmap    []int
}

var _ fem.EssentialBC = &Essential{}

func (e *Essential) AtDof(i int) (bool, float64) {
	dofCount := len(e.dofmap)
	if len(e.dofmap) == 0 {
		e.dofmap, _ = fem.DofMapping(fem.DofAll, e)
		dofCount = len(e.dofmap)
	}
	lN := len(e.Nodes)
	idiv := i / dofCount
	if idiv >= lN {
		return e.Augmented[i-lN*dofCount], 0
	}
	imod := i % len(e.dofmap)
	d := (e.Nodes[idiv] >> e.dofmap[imod]) & 1
	// return (d != 0) != e.AreFree, 0
	if e.AreFree {
		return d == 0, 0
	}
	return d != 0, 0
}

func (e *Essential) NumberOfDofs() int  { return e.D.Count()*len(e.Nodes) + len(e.Augmented) }
func (e *Essential) Dofs() fem.DofsFlag { return e.D }
