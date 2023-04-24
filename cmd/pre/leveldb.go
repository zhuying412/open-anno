package pre

import (
	"encoding/json"
	"fmt"
	"log"
	"open-anno/pkg"
	"regexp"
	"strings"

	"github.com/go-playground/validator/v10"
	"github.com/spf13/cobra"
	"github.com/syndtr/goleveldb/leveldb"
)

type PreLeveldbParam struct {
	Inputs    []string `validate:"required,pathsexists"`
	Outdir    string   `validate:"required"`
	BatchSize int      `validate:"required"`
}

func (this PreLeveldbParam) Valid() error {
	validate := validator.New()
	validate.RegisterValidation("pathsexists", pkg.CheckPathsExists)
	err := validate.Struct(this)
	if err != nil {
		return err
	}
	return nil
}

func (this PreLeveldbParam) Run() error {
	db, err := leveldb.OpenFile(this.Outdir, nil)
	if err != nil {
		return err
	}
	defer db.Close()
	headers := make([]map[string]string, 0)
	batch := new(leveldb.Batch)
	i := 0
	fmt.Println(this.Inputs)
	for _, input := range this.Inputs {
		reader, err := pkg.NewIOReader(input)
		if err != nil {
			return err
		}
		scanner := pkg.NewIOScanner(reader)
		for scanner.Scan() {
			text := scanner.Text()
			if strings.HasPrefix(text, "##INFO") {
				re := regexp.MustCompile(`ID=(.+),Number=(.+),Type=(.+),Description="(.+)"`)
				match := re.FindStringSubmatch(text)
				if match != nil {
					headers = append(headers, map[string]string{"Id": match[1], "Number": match[2], "Type": match[3], "Description": match[4]})
				} else {
					fmt.Println("No match header Info")
				}
			}
			if !strings.HasPrefix(text, "#") {
				row := strings.Split(text, "\t")
				key := fmt.Sprintf("%s:%s:%s:%s", row[0], row[1], row[3], row[4])
				info := make(map[string]string)
				for _, item := range strings.Split(row[7], ";") {
					if strings.Contains(item, "=") {
						val := strings.Split(item, "=")
						info[val[0]] = val[1]
					}
				}
				infoText, err := json.Marshal(info)
				if err != nil {
					return err
				}
				batch.Put([]byte(key), infoText)
				if batch.Len() >= this.BatchSize {
					if err := db.Write(batch, nil); err != nil {
						return err
					}
					batch.Reset()
					fmt.Println(i)
				}
				i++
			}
		}
	}
	if batch.Len() > 0 {
		if err := db.Write(batch, nil); err != nil {
			return err
		}
		print(i)
	}
	headerText, err := json.Marshal(headers)
	if err != nil {
		return err
	}
	if err := db.Put([]byte("header:info"), headerText, nil); err != nil {
		return err
	}
	return nil
}

func NewPreLeveldbCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "leveldb",
		Short: "Prepare LevelDB from VCF",
		Run: func(cmd *cobra.Command, args []string) {
			var param PreLeveldbParam
			param.Inputs, _ = cmd.Flags().GetStringArray("input")
			param.Outdir, _ = cmd.Flags().GetString("outdir")
			param.BatchSize, _ = cmd.Flags().GetInt("batch_size")
			err := param.Valid()
			if err != nil {
				cmd.Help()
				log.Fatal(err)
			}
			err = param.Run()
			if err != nil {
				log.Fatal(err)
			}
		},
	}
	cmd.Flags().StringArrayP("input", "i", []string{}, "Input VCF Files")
	cmd.Flags().StringP("outdir", "o", "", "Output Directory")
	cmd.Flags().IntP("batch_size", "s", 1000*1000, "Parameter Batch Size")
	return cmd
}
