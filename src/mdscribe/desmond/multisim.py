import sys
import io
import os
import errno
import importlib.resources
import pathlib
import pyparsing as pp
from pyparsing import pyparsing_common as ppc
from dotmap import DotMap




class Multisim:
    """Parsing Desmond multisim .cfg and .msj expressions"""
    # variable, value, array, expr
    # pyparsing module’s default behavior is to ignore whitespace.
    # +: AND, |: MatchFirst, left-to-right, ^: Or(longest match)
    # Group: --> list
    # Dict: --> dict
    # Forward: --> recursive

    EQ = pp.Suppress('=')
    LBRACKET, RBRACKET, LBRACE, RBRACE = map(pp.Literal, "[]{}")
    variable = (pp.Word(pp.alphanums + "._/?-@") + 
                pp.Opt("." + pp.Word(pp.alphanums))).set_parse_action(''.join)
    _string1 = pp.Word(pp.alphanums + "._/?-@*")
    _string2 = pp.quoted_string()
    _number  = ppc.number()
    value   = (_string1 | _string2 | _number)
    array   = pp.Forward()
    array   <<= pp.Group(LBRACKET + (pp.ZeroOrMore(value | array)) + RBRACKET)
    expr    = pp.Forward()
    _expr_0  = (variable + EQ + value)
    _expr_1  = (variable + EQ + array)
    _expr_2  = (variable + EQ + pp.Group(
        LBRACE + pp.ZeroOrMore(expr) + RBRACE ))
    _expr_3  = (variable + EQ + pp.Group(
        LBRACKET + pp.ZeroOrMore(pp.Group(LBRACE + expr + RBRACE)) + RBRACKET))
    _expr_4  = pp.Group(variable + pp.Group(
        LBRACE + pp.ZeroOrMore( expr ) + RBRACE))
    expr    <<= pp.OneOrMore(pp.Dict(pp.Group(
        _expr_0 | _expr_1 | _expr_2 | _expr_3 | _expr_4)))
    expr.ignore("#" + pp.restOfLine)


    def __init__(self, **kwargs):
        self.template_path = None
        self.ParseResults = None
        self.dict = {}
        self.dot = DotMap()
        self.output = sys.stdout # do not attempt to close
        self.indent = 4

        if 'template' in kwargs:
            template = kwargs['template']
            template_path = pathlib.Path(template)
            if template_path.is_file():
                self.template_path = template_path
            else:
                with importlib.resources.files('mdscribe.desmond') as template_path:
                    self.template_path = template_path / template
            if self.template_path is None:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), template)

            self.ParseResults = Multisim.expr.parse_file(self.template_path)
            self.decode()

        elif 'string' in kwargs:
            string = kwargs['string']
            try:
                self.ParseResults = Multisim.expr.parse_string(string)
                self.decode()
            except:
                raise RuntimeError("Multisim: cannot parse the input string.")
        else:
            raise RuntimeError("Multisim: template filename or string is required.")



    @staticmethod
    def unfold_dict_key(d) -> None:
        """Unfold '.' in the keys of a dictionary

        Args:
            d (dict): dictionary (to be modified)
        """
        ks = [ k for k in d ]
        for k in ks:
            v = d[k]
            if '.' in k:
                # d[k] = v --> d[kk[0]] = { kk[1] : v }
                _u = None
                _k = k.split('.')[0]
                for kk in k.split('.')[-1:0:-1]:
                    if _u is None:
                        _u = {kk : v }
                    else:
                        _u = {kk : _u}
                d[_k] = _u
                del d[k]
            else:
                pass
            if isinstance(v, dict):
                Multisim.unfold_dict_key(v)


    @staticmethod
    def traverse_dict(d) -> None:
        """Recursively traverse a nested dictionary/list

        Args:
            d (dict): dictionary (to be modified)
        """
        if isinstance(d, dict):
            for k,v in d.items():
                if isinstance(v, dict):
                    Multisim.traverse_dict(v)            
                elif isinstance(v, list):
                    if v == ['{','}']:
                        d[k] = {}
                    elif v == ['[',']']:
                        d[k] = []
                    elif v[0] == '[' and v[-1] == ']':
                        d[k] = v[1:-1]
                        Multisim.traverse_dict(d[k])
                    elif v[0] == '{' and v[-1] == '}':
                        d[k] = dict(v[1:-1])
                        Multisim.traverse_dict(d[k])
        elif isinstance(d, list):
            for k,v in enumerate(d):
                if isinstance(v, dict):
                    Multisim.traverse_dict(v)
                elif isinstance(v, list):
                    if v == ['{','}']:
                        d[k] = {}
                    elif v == ['[',']']:
                        d[k] = []
                    elif v[0] == '[' and v[-1] == ']':
                        d[k] = v[1:-1]
                        Multisim.traverse_dict(d[k])
                    elif v[0] == '{' and v[-1] == '}':
                        d[k] = dict(v[1:-1])
                        Multisim.traverse_dict(d[k])


    def decode(self) -> None:
        """decode the parsed results into a dictionary and its dotmap"""
        # create .dict
        if isinstance(self.ParseResults.as_list()[0][0], str): # key
            self.dict = self.ParseResults.as_dict()
            Multisim.traverse_dict(self.dict)
            # handle the case where key has '.'
            Multisim.unfold_dict_key(self.dict)

        elif isinstance(self.ParseResults.as_list()[0][0], list):
            self.dict = [] # now self.dict is a list of dictionary
            for section in self.ParseResults:
                dict_ = dict(section.as_list())
                Multisim.traverse_dict(dict_)
                # handle the case where key has '.'
                Multisim.unfold_dict_key(dict_)
                self.dict.append(dict_)

        # create .dot
        if isinstance(self.dict, list):
            self.dot = {}
            for d in self.dict:
                for k,v in d.items():
                    if k in self.dot:
                        if isinstance(self.dot[k], list):
                            self.dot[k].append(v)
                        else:
                            self.dot[k] = [self.dot[k], v]
                    else:
                        self.dot[k] = v
            self.dot = DotMap(self.dot)
        else:
            self.dot = DotMap(self.dict)

   

    def write(self, output=None):
        """Writes DOT object"""
        if isinstance(output, io.IOBase):
            self.output = output
        elif isinstance(output, str):
            self.output = open(output, "w")

        if isinstance(self.dict, list):
            blocks = []
            for k, v in self.dot.items():
                if isinstance(v, list):
                    for vv in v:
                        blocks.append({k:vv})
                else:
                    blocks.append({k:v})
            for block in blocks:
                self._write_dict(block, block=True)
                self.output.write("\n")
        else:
            self._write_dict(self.dot)


    def _write_dict(self, d, block=False, depth=0):
        """subroutine of .write() method"""
        spc = ' ' * self.indent * depth
        if isinstance(d, dict) or isinstance(d, DotMap):
            for k, v in d.items():
                k = str(k)
                if v:
                    if isinstance(v, dict) or isinstance(v, DotMap):
                        if depth == 0 and block:
                            self.output.write(spc + k + " {\n")
                        else:
                            self.output.write(spc + k + " = {\n")
                        self._write_dict(v, depth=depth+1)
                        self.output.write(spc + "}\n")
                    elif isinstance(v, list):
                        self.output.write(spc + k + " = [")
                        for vv in v:
                            if isinstance(vv, dict) or isinstance(vv, DotMap): 
                                self.output.write("{\n")
                            elif isinstance(vv, list):
                                self.output.write("[")
                            self._write_dict(vv, depth=depth+1)
                            if isinstance(vv, dict) or isinstance(vv, DotMap):
                                self.output.write(spc + "}")
                            elif isinstance(vv, list):
                                self.output.write("]")                
                        self.output.write("]\n")
                    else:
                        self.output.write(spc + k + " = " + str(v) +"\n")
                else:
                    if isinstance(v, list) and (not bool(v)):
                        self.output.write(spc + k + " = []\n")
                    elif (isinstance(v, dict) or isinstance(v, DotMap))and (not bool(v)):
                        self.output.write(spc + k + " = {\n}\n")
                    else:
                        self.output.write(spc + k + " =   \n")
        elif isinstance(d, list):
            for v in d:
                self._write_dict(v, depth=depth+1)
        else:
            self.output.write(" " + str(d) + " ")


    def to_dot(self) -> DotMap:
        """Returns parsed results as a DotMap object

        Returns:
            DotMap : DotMap object
        """
        return self.dot
    

    def to_list(self) -> list:
        """Returns parsed results as a list

        Returns:
            list : list
        """
        return self.ParseResults.as_list()


    def to_dict(self) -> dict:
        """Returns parsed results as a dictionary

        Returns:
            dict : dictionary
        """
        return self.dict