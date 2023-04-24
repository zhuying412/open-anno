package db

import (
	"encoding/json"
	"errors"
	"open-anno/pkg"

	"github.com/brentp/vcfgo"
	"github.com/syndtr/goleveldb/leveldb"
)

func GetHeaderLevelDB(db *leveldb.DB) ([]*vcfgo.Info, error) {
	val, err := db.Get([]byte("header:info"), nil)
	if err != nil {
		return []*vcfgo.Info{}, err
	}
	var headInfo []*vcfgo.Info
	err = json.Unmarshal(val, &headInfo)
	if err != nil {
		return []*vcfgo.Info{}, err
	}
	return headInfo, nil
}

func AnnoFilterBasedLevelDB(variant pkg.IVariant, db *leveldb.DB) (map[string]any, error) {
	annoInfo := make(map[string]any)
	val, err := db.Get([]byte(variant.PK()), nil)
	if err != nil {
		if errors.Is(err, leveldb.ErrNotFound) {
			return annoInfo, nil
		}
		return map[string]any{}, err
	}
	err = json.Unmarshal(val, &annoInfo)
	if err != nil {
		return map[string]any{}, err
	}
	return annoInfo, nil
}
