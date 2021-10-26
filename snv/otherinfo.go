package snv

import "strings"

type OtherInfo map[string]string

func NewOtherInfo(info string) OtherInfo {
	otherinfo := make(OtherInfo)
	fileds := strings.Split(info, ",")
	for _, filed := range fileds {
		kv := strings.Split(filed, "=")
		if len(kv) > 1 {
			otherinfo[kv[0]] = kv[1]
		}
	}
	return otherinfo
}
