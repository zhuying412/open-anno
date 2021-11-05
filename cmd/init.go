package cmd

import "github.com/spf13/cobra"

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
	RootCmd.AddCommand(NewPrepareCmd())
	RootCmd.AddCommand(NewSnvCmd())
	RootCmd.AddCommand(NewCnvCmd())
	RootCmd.AddCommand(NewIndexDatabaseCmd())
}
