def handle_unknown_scale(scale):
    raise Exception(f"Axis scale '{scale}' unknown.")


def handle_unknown_solver(solver):
    raise Exception(f"Solver '{solver}' unknown.")
