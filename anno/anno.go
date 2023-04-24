package anno

import (
	"open-anno/anno/db"
	"open-anno/anno/gene"
	"open-anno/pkg"
	"sync"

	"github.com/brentp/bix"
	"github.com/brentp/faidx"
	"github.com/syndtr/goleveldb/leveldb"
)

type AnnoInfo struct {
	PK    string
	Data  map[string]any
	Error error
}

func (this *AnnoInfo) AddAnno(anno map[string]any) {
	for key, val := range anno {
		this.Data[key] = val
	}
}

func AnnoSnv(snv *pkg.SNV, gpeTbx *bix.Bix, fbTbxs []*bix.Bix, rbTbxs []*bix.Bix, dbnames []string, fbDBs []*leveldb.DB, genome *faidx.Faidx, overlap float64) AnnoInfo {
	annoInfo := AnnoInfo{PK: snv.PK(), Error: nil, Data: make(map[string]any)}
	var anno map[string]any
	var err error
	// anno, err = gene.AnnoSnv(snv, gpeTbx, genome)
	// if err != nil {
	// 	annoInfo.Error = err
	// 	return annoInfo
	// }
	// annoInfo.AddAnno(anno)
	for _, tbx := range fbTbxs {
		anno, err = db.AnnoFilterBased(snv, tbx)
		if err != nil {
			annoInfo.Error = err
			return annoInfo
		}
		annoInfo.AddAnno(anno)
	}
	for i, tbx := range rbTbxs {
		anno, err = db.AnnoRegionBased(snv, tbx, overlap, dbnames[i])
		if err != nil {
			annoInfo.Error = err
			return annoInfo
		}
		annoInfo.AddAnno(anno)
	}
	for _, leveldb := range fbDBs {
		anno, err = db.AnnoFilterBasedLevelDB(snv, leveldb)
		if err != nil {
			annoInfo.Error = err
			return annoInfo
		}
		annoInfo.AddAnno(anno)
	}
	return annoInfo
}

func AnnoSnvWorker(snvs chan *pkg.SNV, gpeTbx *bix.Bix, fbTbxs []*bix.Bix, rbTbxs []*bix.Bix, dbnames []string, fbDBs []*leveldb.DB, genome *faidx.Faidx, overlap float64, result chan AnnoInfo, wg *sync.WaitGroup) {
	defer wg.Done()
	for snv := range snvs {
		result <- AnnoSnv(snv, gpeTbx, fbTbxs, rbTbxs, dbnames, fbDBs, genome, overlap)
	}
}

func AnnoCnv(cnv *pkg.CNV, gpeTbx *bix.Bix, gbTbxs []*bix.Bix, dbnames []string, overlap float64) AnnoInfo {
	annoInfo := AnnoInfo{PK: cnv.PK(), Error: nil, Data: make(map[string]any)}
	anno, err := gene.AnnoCnv(cnv, gpeTbx)
	if err != nil {
		annoInfo.Error = err
		return annoInfo
	}
	annoInfo.AddAnno(anno)
	for i, tbx := range gbTbxs {
		anno, err = db.AnnoRegionBased(cnv, tbx, overlap, dbnames[i])
		if err != nil {
			annoInfo.Error = err
			return annoInfo
		}
		annoInfo.AddAnno(anno)
	}
	return annoInfo
}

func AnnoCnvWorker(cnvs chan *pkg.CNV, gpeTbx *bix.Bix, rbTbxs []*bix.Bix, dbnames []string, overlap float64, result chan AnnoInfo, wg *sync.WaitGroup) {
	defer wg.Done()
	for cnv := range cnvs {
		result <- AnnoCnv(cnv, gpeTbx, rbTbxs, dbnames, overlap)
	}
}
