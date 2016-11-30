typedef char[2] twoByte;


int main() {
	static_assert(2*sizeof(char) == sizeof(twoByte), "diff sizes");
}
