from reservoirflow.base import Base

class Plot(Base):
    
    name = "Plot"
    
    def __init__(self):
        pass
        
        
if __name__ == "__main__":
    plot = Plot()
    print(plot)