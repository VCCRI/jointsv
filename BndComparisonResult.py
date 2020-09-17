class BndComparisonResult:
    """
    This object represents the result of a SV call operation
    """
    def __init__(self, is_sv, svtype, initial_position, final_position):
        self.is_sv = is_sv
        self.svtype = svtype
        self.initial_position = initial_position
        self.final_position = final_position
