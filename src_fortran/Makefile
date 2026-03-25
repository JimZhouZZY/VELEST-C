FC = gfortran
FFLAGS = -std=legacy
TARGET = velest

# Expand the wildcard to get actual files
SOURCES = velest.f $(wildcard sub/*.f)

$(TARGET): $(SOURCES)
	$(FC) $(FFLAGS) $^ -o $@

clean:
	rm -f $(TARGET)

.PHONY: clean
