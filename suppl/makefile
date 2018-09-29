all: B
B: B.c
	gcc B.c -lm -o B
clean:
	rm -f ./B ./B.h

# a not very elegant way to remove conditions with empty bodies
EMPTYBODYSED=sed -e 's/\r[^\r]*{\r[ ]*}//g'
EMPTYBODY=tr '\n' '\r' | $(EMPTYBODYSED) | $(EMPTYBODYSED) | tr '\r' '\n'
B.h : B
	./B f  | $(EMPTYBODY) > B.h
