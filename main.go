package main

import (
	"log"
	"open-anno/cmd"
	"open-anno/cmd/anno"
	"open-anno/cmd/pre"
	"open-anno/cmd/tools"

	"github.com/spf13/cobra"
)

var RootCmd *cobra.Command

func NewRootCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "OpenAnno",
		Short: "A Genome Annotate Tool",
	}
	return cmd
}

func NewPreCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "pre",
		Short: "Prepare database",
	}
	cmd.AddCommand(pre.NewPreGeneBasedCmd())
	cmd.AddCommand(pre.NewIndexDatabaseCmd())
	cmd.AddCommand(pre.NewGeneInfoCmd())
	return cmd
}

func NewAnnoCmd() *cobra.Command {

	annoSnv := &cobra.Command{
		Use:   "snv",
		Short: "Annotate SNV",
	}
	annoSnv.AddCommand(anno.NewAnnoCmd("gb", anno.VType_SNV, anno.DType_G))
	annoSnv.AddCommand(anno.NewAnnoCmd("fb", anno.VType_SNV, anno.DType_F))
	annoSnv.AddCommand(anno.NewAnnoCmd("rb", anno.VType_SNV, anno.DType_R))
	annoCnv := &cobra.Command{
		Use:   "cnv",
		Short: "Annotate SNV",
	}
	annoCnv.AddCommand(anno.NewAnnoCmd("gb", anno.VType_CNV, anno.DType_G))
	annoCnv.AddCommand(anno.NewAnnoCmd("rb", anno.VType_CNV, anno.DType_R))
	cmd := &cobra.Command{
		Use:   "anno",
		Short: "Annotate variants",
	}
	cmd.AddCommand(annoSnv)
	cmd.AddCommand(annoCnv)
	cmd.AddCommand(anno.NewMergeCmd())
	return cmd
}

func NewToolsCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "tools",
		Short: "Tools",
	}
	cmd.AddCommand(tools.NewVCF2AnnoInputCmd())
	cmd.AddCommand(tools.NewAnnoInput2VCFCmd())
	return cmd
}

func init() {
	RootCmd = NewRootCmd()
	RootCmd.AddCommand(NewPreCmd())
	RootCmd.AddCommand(NewAnnoCmd())
	RootCmd.AddCommand(NewToolsCmd())
	RootCmd.AddCommand(cmd.NewTestCmd())
}

func main() {
	err := RootCmd.Execute()
	if err != nil {
		log.Panic(err)
	}
}
