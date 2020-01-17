package core

func IfElse(condition bool, ifValue interface{}, elseValue interface{}) interface{} {
	if condition {
		return ifValue
	}
	return elseValue
}
