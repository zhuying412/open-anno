package command

import (
	"OpenAnno/command/anno"
	"OpenAnno/command/prepare"
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

func init() {
	RootCmd = NewRootCmd()
	RootCmd.AddCommand(prepare.NewPrepareCmd())
	RootCmd.AddCommand(anno.NewAnnoCmd())
}
