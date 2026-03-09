using System.IO;
using Sharpmake;

namespace Qaws
{
    [Generate]
    public class QawsLibProject : CommonProject
    {
        public QawsLibProject()
        {
            Name = "Qaws";
            SourceRootPath = @"[project.SharpmakeCsPath]\..\..\..\src\qaws";
        }

        [Configure()]
        public override void ConfigureAll(Configuration conf, QawsTarget target)
        {
            base.ConfigureAll(conf, target);

            // Output type: Static library
            conf.Output = Configuration.OutputType.Lib;

            // Include paths
            conf.IncludePaths.Add(@"[project.SharpmakeCsPath]\..\..\..\src\qaws");

            // Source files (explicitly list them for now)
            conf.SourceFilesBuildExcludeRegex.Add(@".*\.sharpmake\.cs$");
        }
    }
}
