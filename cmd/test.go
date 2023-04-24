package cmd

import (
	"fmt"
	"log"
	"sync"
	"time"

	"github.com/spf13/cobra"
)

func Test(in chan int, out chan int, wg *sync.WaitGroup) {
	defer wg.Done()
	for i := range in {
		fmt.Printf("%d start\n", i)
		time.Sleep(3 * time.Second)
		fmt.Printf("%d end\n", i)
		out <- i * i
	}

}

func RunTest() {
	in := make(chan int, 100)
	out := make(chan int, 100)
	// 将要处理的元素发送到通道中
	for i := 1; i <= 100; i++ {
		in <- i
	}
	close(in)
	var wg sync.WaitGroup
	for i := 0; i < 3; i++ {
		wg.Add(1)
		go Test(in, out, &wg)
		// go func() {
		// 	for {
		// 		if _, ok := <-in; !ok {
		// 			break
		// 		}
		// 		Test(in, out, &wg)

		// 	}
		// }()
	}
	// wg.Wait()
	// close(ch)
	go func() {
		// 阻塞等待所有的 worker 结束
		wg.Wait()
		// 所有 worker 结束后，关闭 resultChan 通道
		close(out)
	}()
	// close(ch)
	for i := range out {
		fmt.Println(i)
	}
}

func NewTestCmd() *cobra.Command {
	cmd := &cobra.Command{
		Use:   "test",
		Short: "Test",
		Run: func(cmd *cobra.Command, args []string) {
			input, _ := cmd.Flags().GetString("input")
			if input == "" {
				err := cmd.Help()
				if err != nil {
					log.Fatal(err)
				}
			} else {
				RunTest()
			}
		},
	}
	cmd.Flags().StringP("input", "i", "", "input")
	return cmd
}
