/*  Copyright (C) 2004 Glenn Pierce.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "itkFDFCommonImageIO.h"

// Remove a particular type of character from a string
void RemoveCharacters(std::string &line, std::string character)
{
    std::string::size_type characterPosition = line.find_first_of(character, 0);
                                                                                                                                while ( characterPosition != std::string::npos ) {
        line.erase(characterPosition, 1);
        characterPosition = line.find_first_of(character, characterPosition);
    }
}

void Tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters)
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

std::string ParseLine(std::string line)
{
    // strip *
    RemoveCharacters(line, "*");
    RemoveCharacters(line, "\"");
    RemoveCharacters(line, "[");
    RemoveCharacters(line, "]");

    // Need to deal with space between {}
    std::string::size_type startBracketPosition = line.find_first_of("{", 0);
    std::string::size_type endBracketPosition = line.find_first_of("}", startBracketPosition);
                                                                                                                            
    if ( startBracketPosition != std::string::npos && endBracketPosition != std::string::npos) {
        std::string element = line.substr(startBracketPosition, endBracketPosition - startBracketPosition);
                                                                                                                            
        // Find whitespace within {} and erase
        std::string::size_type whiteSpacePosition = line.find_first_of(" ", startBracketPosition);
                                                                                                                            
        while (whiteSpacePosition != std::string::npos)
        {
            line.erase(whiteSpacePosition, 1);
            whiteSpacePosition = line.find_first_of(" ", whiteSpacePosition);
        }
    }

    return line;
}
