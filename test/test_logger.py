import ROOT
import supera


def test_logger_uses_configured_default_level():
    original_level = ROOT.supera.Logger.default_level()
    try:
        ROOT.supera.Logger.default_level(ROOT.supera.msg.kWARNING)

        direct_logger = ROOT.supera.Logger("direct")
        named_logger = ROOT.supera.Logger.get("named")

        assert direct_logger.level() == ROOT.supera.msg.kWARNING
        assert named_logger.level() == ROOT.supera.msg.kWARNING
    finally:
        ROOT.supera.Logger.default_level(original_level)
