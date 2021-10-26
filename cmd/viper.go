package cmd

import (
	"github.com/spf13/viper"
	"log"
)

func InitViper(dbPath string) {
	viper.AddConfigPath(dbPath)
	viper.SetConfigName("config")
	err := viper.ReadInConfig() // 根据以上配置读取加载配置文件
	if err != nil {
		log.Fatal(err) // 读取配置文件失败致命错误
	}
	viper.Set("db.path", dbPath)
}
