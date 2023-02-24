package variant

import (
	"fmt"
)

const (
	VType_SNP = "SNP"
	VType_INS = "INS"
	VType_DEL = "DEL"
	VType_DUP = "DUP"
	VType_SUB = "SUB"
)

type AnnoVariant struct {
	Chrom string `json:"chrom"`
	Start int    `json:"start"`
	End   int    `json:"end"`
	Ref   string `json:"ref"`
	Alt   string `json:"alt"`
}

func (this AnnoVariant) Type() string {
	if this.Ref == "DIP" {
		return this.Alt
	}
	if this.Ref == "-" {
		return VType_INS
	} else if this.Alt == "-" {
		return VType_DEL
	} else {
		if len(this.Ref) > 1 || len(this.Alt) > 1 {
			return VType_SUB
		}
		return VType_SNP
	}
}

func (this AnnoVariant) PK() string {
	return fmt.Sprintf("%s:%d:%d:%s:%s", this.Chrom, this.Start, this.End, this.Ref, this.Alt)
}

func (this AnnoVariant) Less(that AnnoVariant) bool {
	if this.Chrom < that.Chrom {
		return true
	}
	if this.Chrom == that.Chrom {
		if this.Start < that.Start {
			return true
		}
		if this.Start == that.Start {
			if this.End < that.End {
				return true
			}
			if this.End == that.End {
				if this.Ref < that.Ref {
					return true
				}
				if this.Ref == that.Ref {
					if this.Alt < that.Alt {
						return true
					}
				}
			}
		}
	}
	return false
}

func (this AnnoVariant) Equal(that AnnoVariant) bool {
	return this.Chrom == that.Chrom && this.Start == that.Start && this.End == that.End && this.Ref == that.Ref && this.Alt == that.Alt
}

func (this AnnoVariant) Greater(that AnnoVariant) bool {
	if this.Chrom > that.Chrom {
		return true
	}
	if this.Chrom == that.Chrom {
		if this.Start > that.Start {
			return true
		}
		if this.Start == that.Start {
			if this.End > that.End {
				return true
			}
			if this.End == that.End {
				if this.Ref > that.Ref {
					return true
				}
				if this.Ref == that.Ref {
					if this.Alt > that.Alt {
						return true
					}
				}
			}
		}
	}
	return false
}

type IVariant interface {
	AnnoVariant() AnnoVariant
}
