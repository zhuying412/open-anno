package main

import (
	"github.com/spf13/cobra"
	"grandanno/cnv"
	"grandanno/db"
	"grandanno/snv"
	"log"
)

var rootCmd *cobra.Command

func NewRootCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "OpenAnno",
		Short: "A Genome Annotate Tool",
	}
	return cmd
}

func NewPrepareCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "prepare",
		Short: "Prepare required database files",
		Run: func(cmd *cobra.Command, args []string) {
			log.Println(cmd.Flag("database").Value.String())
			if cmd.Flag("database").Value.String() == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err.Error())
				}
			} else {
				db.PrepareDatabase(cmd.Flag("database").Value.String(), cmd.Flag("genome").Value.String())
			}
		},
	}
	cmd.Flags().StringP("database", "d", "", "Database Path")
	cmd.Flags().StringP("genome", "g", "hg19", "Genome Version")
	return cmd
}

func NewSnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "snv",
		Short: "Annotate SNV",
		Run: func(cmd *cobra.Command, args []string) {
			log.Println(cmd.Flag("database").Value.String())
			if cmd.Flag("database").Value.String() == "" ||
				cmd.Flag("input").Value.String() == "" ||
				cmd.Flag("output").Value.String() == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err.Error())
				}
			} else {
				snv.AnnoSnv(
					cmd.Flag("database").Value.String(),
					cmd.Flag("genome").Value.String(),
					cmd.Flag("input").Value.String(),
					cmd.Flag("output").Value.String(),
				)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "SNV input file")
	cmd.Flags().StringP("output", "o", "", "Annotated ouput file")
	cmd.Flags().StringP("database", "d", "", "Database Path")
	cmd.Flags().StringP("genome", "g", "hg19", "Genome Version")
	return cmd
}

func NewCnvCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "cnv",
		Short: "Annotate CNV",
		Run: func(cmd *cobra.Command, args []string) {
			log.Println(cmd.Flag("database").Value.String())
			if cmd.Flag("database").Value.String() == "" ||
				cmd.Flag("input").Value.String() == "" ||
				cmd.Flag("output").Value.String() == "" {
				err := cmd.Help()
				if err != nil {
					log.Panic(err.Error())
				}
			} else {
				cnv.AnnoCnv(
					cmd.Flag("database").Value.String(),
					cmd.Flag("genome").Value.String(),
					cmd.Flag("input").Value.String(),
					cmd.Flag("output").Value.String(),
				)
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "SNV input file")
	cmd.Flags().StringP("output", "o", "", "Annotated ouput file")
	cmd.Flags().StringP("database", "d", "", "Database Path")
	cmd.Flags().StringP("genome", "g", "hg19", "Genome Version")
	return cmd
}

func init() {
	rootCmd = NewRootCmd()
	rootCmd.AddCommand(NewPrepareCmd())
	rootCmd.AddCommand(NewSnvCmd())
	rootCmd.AddCommand(NewCnvCmd())
}

func main() {
	err := rootCmd.Execute()
	if err != nil {
		log.Panic(err.Error())
	}
}
