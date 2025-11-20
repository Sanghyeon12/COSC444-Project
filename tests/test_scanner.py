# test_scanner.py - Sadman

def test_scan_sequence_empty():
    from src.scan.scanner import scan_sequence
    assert scan_sequence("id", "", {}) == []
